/*
 * Copyright (c) 1998 - 2010. University Corporation for Atmospheric
 * Research/Unidata
 * Portions of this software were developed by the Unidata Program at the
 * University Corporation for Atmospheric Research.
 *
 * Access and use of this software shall impose the following obligations
 * and understandings on the user. The user is granted the right, without
 * any fee or cost, to use, copy, modify, alter, enhance and distribute
 * this software, and any derivative works thereof, and its supporting
 * documentation for any purpose whatsoever, provided that this entire
 * notice appears in all copies of the software, derivative works and
 * supporting documentation.  Further, UCAR requests that the user credit
 * UCAR/Unidata in any publications that result from the use of this
 * software or in any product that includes this software. The names UCAR
 * and/or Unidata, however, may not be used in any advertising or publicity
 * to endorse or promote any products or commercial entity unless specific
 * written permission is obtained from UCAR/Unidata. The user also
 * understands that UCAR/Unidata is not obligated to provide the user with
 * any support, consulting, training or assistance of any kind with regard
 * to the use, operation and performance of this software nor to provide
 * the user with any updates, revisions, new versions or "bug fixes."
 *
 * THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
 * INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
 * FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
 * WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
 */


package ucar.nc2.iosp.sigmet;

//~--- non-JDK imports --------------------------------------------------------

import ucar.nc2.time.CalendarDate;
import ucar.unidata.io.RandomAccessFile;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

/**
 * @author yuanho
 * @since Apr 7, 2010
 */
public class SigmetVolumeScan {
    private static String[] dataTypeNames = {"Extended Headers", "TotalPower",
            "Reflectivity", "Velocity", "Width", "DifferentialReflectivity",
            "CorrectedReflectivity", "TotalPower", "Reflectivity", "Velocity",
            "Width", "DifferentialReflectivity", "RainfallRate",
            "SpecificDifferentialPhase", "SpecificDifferentialPhase",
            "DifferentialPhase", "CorrectedVelocity", "SQI",
            "CorrelationCoefficient", "CorrelationCoefficient",
            "CorrectedReflectivity", "CorrectedVelocity", "SQI",
            "DifferentialPhase", "LinearDepolarizationRatioH",
            "LinearDepolarizationRatioH", "LinearDepolarizationRatioV",
            "LinearDepolarizationRatioV"};

    private static String[] dataTypeUnit = {"", "dBZ", "dBZ", "m/sec",
            "m/sec", "dB", "dBZ", "dBZ", "dBZ", "m/sec", "m/sec", "dB", "mm/hr",
            "deg/km", "deg/km", "deg", "m/sec", "dimensionless",
            "dimensionless", "dimensionless", "dBZ", "m/sec", "dimensionless",
            "deg", "dB", "dB", "dB", "dB"};

    private Map<String, List<List<Ray>>> groups;
    private Map<Integer, List<Ray>> buffers;
    private int[] num_gates;
    public CalendarDate[] date;
    CalendarDate start_date, end_date;
    public Ray firstRay = null;
    public Ray lastRay = null;
    public ucar.unidata.io.RandomAccessFile raf;

    static final int REC_SIZE = 6144;
    private RecordsHeader recHdr;

    public class RecordsHeader {
        /**
         * Read some global data from SIGMET file. The SIGMET file consists of
         * records with
         * fixed length=6144 bytes.
         */
        private int prf, wave;
        public float vNyq, velScale;

        public short ground_height, radar_height;
        public float radar_lat, radar_lon;
        public int radar_alt;

        public short num_rays, multiprf, num_bins, num_sweeps;
        public float range_first, range_last;
        public int step;

        public CalendarDate base_date;

        public String stnName, stnName_util;
        public List<Integer> paramIds;
        public List<String> params;
        public List<String> paramUnits;

        public RecordsHeader(ucar.unidata.io.RandomAccessFile raf) {
            try {
                //      -- Read from <product_end> of the 1st record -- 12+320+120
                //      -- Calculate Nyquist velocity --------------------
                raf.seek(452);
                prf = raf.readInt();
                raf.seek(480);
                wave = raf.readInt();
                vNyq = calcNyquist(prf, wave);

                //      -- Read from the 2nd record----------- 6144+12
                // (struct_hdr)+168(from ingest_config)
                raf.seek(6288);
                stnName = raf.readString(16);
                raf.seek(6306);
                stnName_util = raf.readString(16);

                raf.skipBytes(4);
//                raf.seek(6324);
                radar_lat = calcAngle(raf.readInt());
                radar_lon = calcAngle(raf.readInt()); // 6328
                ground_height = raf.readShort(); // 6332
                radar_height = raf.readShort(); // 6334

                raf.skipBytes(4);
                num_rays = raf.readShort(); // 6340

                raf.skipBytes(2);
                radar_alt = raf.readInt(); // 6344

                // Determine what moments are available in the file
                raf.seek(6772);
                int data_mask = raf.readInt();
                paramIds = new ArrayList<>();
                params = new ArrayList<>();
                paramUnits = new ArrayList<>();
                for (int param = 0; param < 32; ++param) {
                    if (((data_mask >> param) & 0x1) != 0) {
                        paramIds.add(param);
                        params.add(dataTypeNames[param]);
                        paramUnits.add(dataTypeUnit[param]);
                    }
                }

                raf.seek(6912);
                multiprf = raf.readShort();
                velScale = vNyq * (multiprf + 1);

                raf.seek(7408);
                range_first = raf.readInt() * 0.01f; //  cm    7408
                range_last = raf.readInt() * 0.01f; //  cm  7412

                raf.skipBytes(2);
                num_bins = raf.readShort();             //7418
                if (num_bins % 2 != 0)
                    num_bins = (short) (num_bins + 1);

                raf.skipBytes(4);
                step = raf.readInt(); //  cm    7424

                raf.seek(7574);
                num_sweeps = raf.readShort(); // 7574

                raf.seek(12312);
                base_date = readDate(raf);
            } catch (Exception e) {
                System.out.println(e.toString());
                e.printStackTrace();
            }
        }
    }

    public RecordsHeader header() { return recHdr; }

    /**
     * Read all the values from SIGMET-IRIS file which are necessary to fill in
     * the ncfile.
     *
     * @param raf    ucar.unidata.io.RandomAccessFile corresponds to SIGMET
     *               datafile.
     */
    SigmetVolumeScan(ucar.unidata.io.RandomAccessFile raf)
            throws java.io.IOException {
        int len = 2 * REC_SIZE;    // ---- Read from the 3d record-----------
        short nrec = 0,
                nsweep = 1;
        int nwords = 0,
                end_words = 0,
                data_read = 0,
                num_zero = 0,
                rays_count = 0,
                nb = 0,
                pos = 0,
                pos_ray_hdr = 0;
        short a0 = 0,
                a00 = 0;
        short beg_az = 0,
                end_az = 0,
                end_elev = 0,
                num_bins = 0,
                time_start_sw = 0;
        float az, elev, step;
        //     byte      data          = 0;
        boolean beg_rec = true,
                end_rec = true,
                read_ray_hdr = true;
        int cur_len,
                beg = 1,
                kk = 0,
                nu = 0,
                dty;

        // Input
        this.raf = raf;
        raf.order(RandomAccessFile.LITTLE_ENDIAN);

        int fileLength = (int) raf.length();
        recHdr = new RecordsHeader(raf);
        start_date = recHdr.base_date;

        int nparams = recHdr.params.size();
        short number_sweeps = recHdr.num_sweeps;

        float range_first = recHdr.range_first;
        float range_last = recHdr.range_last;

        short[] num_sweep = new short[nparams];
        short[] num_rays_swp = new short[nparams];
        short[] indx_1ray = new short[nparams];
        short[] num_rays_act = new short[nparams];
        short[] angl_swp = new short[nparams];
        short[] bin_len = new short[nparams];
        short[] data_type = new short[nparams];
        num_gates = new int[number_sweeps];
        date = new CalendarDate[nparams * number_sweeps];
        // Array of Ray objects is 2D. Number of columns=number of rays
        // Number of raws = number of types of data if number_sweeps=1,
        // or number of raws = number_sweeps
        groups = new HashMap<>();
        buffers = new HashMap<>();
        List<Ray> time = new ArrayList<>();
        Ray ray = null;

        while (len < fileLength) {
            int rayoffset = 0;
            int rayoffset1 = 0;
            int datalen = 0;

            cur_len = len;

            if (nsweep == number_sweeps & rays_count == beg) {
                return;
            }

            if (beg_rec) {

                // --- <raw_prod_bhdr>  12bytes -----------
                raf.seek(cur_len);
                nrec = raf.readShort();    // cur_len
                nsweep = raf.readShort();    // cur_len+2

                // ---- end of <raw_prod_bhdr> -------------
                cur_len = cur_len + 12;
                beg_rec = false;
            }

            if ((nsweep <= number_sweeps) & (rays_count % beg == 0)) {

                // --Read <ingest_data_hdrs> Number of
                // them=nparams*number_sweeps -----
                // ---Len of <ingest_data_hdr>=76 bytes -----------------
                beg = 0;

                for (int i = 0; i < nparams; i++) {
                    int idh_len = cur_len + 12 + i * 76;

                    raf.seek(idh_len);

                    // Read seconds since midnight
                    date[nu] = readDate(raf);
                    nu++;
                    num_sweep[i] = raf.readShort();     // idh_len+12
                    num_rays_swp[i] = raf.readShort();     // idh_len+14
                    indx_1ray[i] = raf.readShort();     // idh_len+16
                    raf.skipBytes(2);
                    num_rays_act[i] = raf.readShort();
                    beg += num_rays_act[i];    // idh_len+20
                    angl_swp[i] = raf.readShort();     // idh_len+22
                    bin_len[i] = raf.readShort();     // idh_len+24
                    data_type[i] = raf.readShort();     // idh_len+26
                }
                end_date = date[nu - 1];
                cur_len = cur_len + nparams * 76;
            }

            if (end_rec) {

                // --- Read compression code=2 bytes from cur_len
                raf.seek(cur_len);
                a0 = raf.readShort();
                cur_len = cur_len + 2;

                // --- Check if the code=1 ("1" means an end of a ray)
                if (a0 == (short) 1) {
                    if (cur_len % REC_SIZE == 0) {
                        beg_rec = true;
                        end_rec = true;
                        rays_count++;
                        read_ray_hdr = true;
                        pos = 0;
                        data_read = 0;
                        nb = 0;
                        len = cur_len;
                    }
                    else {
                        end_rec = true;
                        len = cur_len;
                        rays_count++;
                    }

                    continue;
                }

                nwords = a0 & 0x7fff;
                end_words = nwords - 6;
                data_read = end_words * 2;
                end_rec = false;

                if (cur_len % REC_SIZE == 0) {
                    len = cur_len;
                    read_ray_hdr = true;
                    beg_rec = true;

                    continue;
                }
            }

            len = cur_len;

            // ---Define output data files for each data_type (= nparams)
            // /sweep ---------
            dty = data_type[0];

            if (nparams > 1) {
                kk = rays_count % nparams;
                dty = data_type[kk];
            }
            else if (number_sweeps > 1) {
                kk = nsweep - 1;
            }

            // --- read ray_header (size=12 bytes=6 words)
            // ---------------------------------------
            if (read_ray_hdr) {
                if (pos_ray_hdr < 2) {
                    raf.seek(cur_len);
                    beg_az = raf.readShort();
                    cur_len = cur_len + 2;
                    len = cur_len;

                    if (cur_len % REC_SIZE == 0) {
                        pos_ray_hdr = 2;
                        beg_rec = true;
                        read_ray_hdr = true;

                        continue;
                    }
                }
                if (pos_ray_hdr < 4) {
                    raf.seek(cur_len);
                    cur_len = cur_len + 2;
                    len = cur_len;

                    if (cur_len % REC_SIZE == 0) {
                        pos_ray_hdr = 4;
                        beg_rec = true;
                        read_ray_hdr = true;
                        continue;
                    }
                }
                if (pos_ray_hdr < 6) {
                    raf.seek(cur_len);
                    end_az = raf.readShort();
                    cur_len = cur_len + 2;
                    len = cur_len;

                    if (cur_len % REC_SIZE == 0) {
                        pos_ray_hdr = 6;
                        beg_rec = true;
                        read_ray_hdr = true;
                        continue;
                    }
                }
                if (pos_ray_hdr < 8) {
                    raf.seek(cur_len);
                    end_elev = raf.readShort();
                    cur_len = cur_len + 2;
                    len = cur_len;

                    if (cur_len % REC_SIZE == 0) {
                        pos_ray_hdr = 8;
                        beg_rec = true;
                        read_ray_hdr = true;

                        continue;
                    }
                }

                if (pos_ray_hdr < 10) {
                    raf.seek(cur_len);
                    num_bins = raf.readShort();
                    cur_len = cur_len + 2;
                    len = cur_len;

                    if (num_bins % 2 != 0) {
                        num_bins = (short) (num_bins + 1);
                    }
                    num_gates[nsweep - 1] = (int) num_bins;
                    if (cur_len % REC_SIZE == 0) {
                        pos_ray_hdr = 10;
                        beg_rec = true;
                        read_ray_hdr = true;

                        continue;
                    }
                }

                if (pos_ray_hdr < 12) {
                    raf.seek(cur_len);
                    time_start_sw = raf.readShort();
                    cur_len = cur_len + 2;
                    len = cur_len;
                }
            }

            // ---------- end of ray header
            // ----------------------------------------------
            az = calcAz(beg_az, end_az);
            elev = calcAngle(end_elev);
            step = SigmetIOServiceProvider.calcStep(range_first, range_last,
                    num_bins);

            if (cur_len % REC_SIZE == 0) {
                len = cur_len;
                beg_rec = true;
                read_ray_hdr = false;

                continue;
            }

            if (pos > 0) {
                data_read = data_read - pos;
                pos = 0;
            }
            if (data_read > 0) {
                raf.seek(cur_len);
                rayoffset = cur_len;
                datalen = data_read;

                for (int i = 0; i < data_read; i++) {
                    cur_len++;
                    nb++;

                    if (cur_len % REC_SIZE == 0) {
                        pos = i + 1;
                        beg_rec = true;
                        read_ray_hdr = false;
                        len = cur_len;
                        raf.seek(cur_len);
                        break;
                    }
                }
                raf.seek(cur_len);
                if (pos > 0) {
                    continue;
                }
            }

            if (cur_len % REC_SIZE == 0) {
                pos = 0;
                beg_rec = true;
                read_ray_hdr = false;
                data_read = 0;
                len = cur_len;

                continue;
            }

            raf.seek(cur_len);
            rayoffset1 = cur_len;

            while (nb < (int) num_bins) {
                a00 = raf.readShort();
                cur_len = cur_len + 2;

                // --- Check if the code=1 ("1" means an end of a ray)
                if (a00 == (short) 1) {
                    // for (int uk = 0; uk < (int) num_bins; uk++) {
                    //  dd[uk] = -999.99f;
                    //  }
                    ray = new Ray(-999.99f, -999.99f, -999.99f, -999.99f,
                            num_bins, (short) (-99), -999, 0, -999,
                            nsweep, dty);
                    rays_count++;
                    beg_rec = false;
                    end_rec = true;

                    break;
                }

                if (a00 < 0) {    // -- This is data
                    nwords = a00 & 0x7fff;
                    data_read = nwords * 2;

                    if (cur_len % REC_SIZE == 0) {
                        pos = 0;
                        beg_rec = true;
                        end_rec = false;
                        len = cur_len;
                        read_ray_hdr = false;

                        break;
                    }

                    raf.seek(cur_len);

                    for (int ii = 0; ii < data_read; ii++) {
                        //   data    = raf.readByte();
                        //    dd[nb]  = SigmetIOServiceProvider.calcData
                        // (recHdr, dty, data);
                        cur_len = cur_len + 1;
                        nb = nb + 1;

                        if (cur_len % REC_SIZE == 0) {
                            pos = ii + 1;
                            beg_rec = true;
                            end_rec = false;
                            len = cur_len;
                            read_ray_hdr = false;
                            raf.seek(cur_len);
                            break;
                        }
                    }
                    raf.seek(cur_len);
                    if (pos > 0) {
                        break;
                    }
                }
                else if (a00 > 0 & a00 != 1) {
                    num_zero = a00 * 2;

                    nb = nb + num_zero;

                    if (cur_len % REC_SIZE == 0) {
                        beg_rec = true;
                        end_rec = false;
                        read_ray_hdr = false;
                        pos = 0;
                        data_read = 0;
                        len = cur_len;

                        break;
                    }
                }
            }                     // ------ end of while for
            // num_bins---------------------------------

            if (cur_len % REC_SIZE == 0) {
                len = cur_len;

                continue;
            }

            raf.seek(cur_len);

            if (nb == (int) num_bins) {
                a00 = raf.readShort();
                cur_len = cur_len + 2;
                end_rec = true;
                ray = new Ray(range_first, step, az, elev, num_bins,
                        time_start_sw, rayoffset, datalen, rayoffset1,
                        nsweep, dty);
                rays_count++;
                if ((nsweep == number_sweeps) & (rays_count % beg == 0)) {
                    addRay(dty, ray);
                    break;
                }

                if (cur_len % REC_SIZE == 0) {
                    beg_rec = true;
                    end_rec = true;
                    read_ray_hdr = true;
                    pos = 0;
                    data_read = 0;
                    nb = 0;
                    len = cur_len;
                    addRay(dty, ray);
                    continue;
                }
            }

            if (firstRay == null) firstRay = ray;

            addRay(dty, ray);

            pos = 0;
            data_read = 0;
            nb = 0;
            read_ray_hdr = true;
            pos_ray_hdr = 0;

            if ((nsweep <= number_sweeps) & (rays_count % beg == 0)) {
                beg_rec = true;
                end_rec = true;
                rays_count = 0;
                nb = 0;
                cur_len = REC_SIZE * (nrec + 1);
                len = cur_len;
                read_ray_hdr = true;
            }

            len = cur_len;
        }    // ------------end of outer while  ---------------
        lastRay = ray;

        for (Integer dataType : buffers.keySet()) {
            groups.put(dataTypeNames[dataType],
                    sortScans(buffers.get(dataType), 1000));
        }
        groups.put("Times", sortScans(time, 1000));
    }

    private void addRay(int dataType, Ray ray)
    {
        if (!buffers.containsKey(dataType))
            buffers.put(dataType, new ArrayList<Ray>());
        buffers.get(dataType).add(ray);
    }
        // --------- fill all of values in the ncfile ------
        // ----------- end of doData -----------------------

    private int max_radials = 0;
    private int min_radials = Integer.MAX_VALUE;

    private List<List<Ray>> sortScans(List<Ray> scans, int siz) {

        // now group by elevation_num
        Map<Short, List<Ray>> groupHash = new HashMap<>(siz);

        for (Ray ray : scans) {
            List<Ray> group = groupHash.get((short) ray.nsweep);

            if (null == group) {
                group = new ArrayList<>();
                groupHash.put((short) ray.nsweep, group);
            }

            group.add(ray);
        }

        for (Short aShort : groupHash.keySet()) {
            List<Ray> group = groupHash.get(aShort);
            Ray[] rr = new Ray[group.size()];
            group.toArray(rr);
            checkSort(rr);
        }

        // sort the groups by elevation_num
        List<List<Ray>> groups = new ArrayList<>(groupHash.values());
        Collections.sort(groups, new GroupComparator());

        // use the maximum radials
        for (List<Ray> group : groups) {
            max_radials = Math.max(max_radials, group.size());
            min_radials = Math.min(min_radials, group.size());
        }

        return groups;
    }

    private class GroupComparator implements Comparator<List<Ray>> {
        public int compare(List<Ray> group1, List<Ray> group2) {
            Ray record1 = group1.get(0);
            Ray record2 = group2.get(0);

            // if (record1.elevation_num != record2.elevation_num)
            return record1.nsweep - record2.nsweep;

            // return record1.cut - record2.cut;
        }
    }

    public List<List<Ray>> getData(String dataType) {
        return groups.get(dataType);
    }

    public int[] getNumberGates() {
        return num_gates;
    }

    public CalendarDate[] getStartSweep() {
        return date;
    }

    /**
     * Sort Ray objects in the same sweep according to the ascended azimuth
     * (from 0 to 360)
     * and time.
     *
     * @param r the array of Ray objects in a sweep. Its length=number_rays
     */
    void checkSort(Ray[] r) {
        int j = 0, n = 0, n1 = 0, n2 = 0;
        short time1 = 0, time2 = 0;
        int[] k1 = new int[300];
        int[] k2 = new int[300];
        //      define the groups of rays with the same "time". For ex.:
        //      group1 - ray[0]={time=1,az=344}, ray[1]={time=1,az=345},
        // ... ray[11]={time=1,az=359}
        //      group2 - ray[12]={time=1,az=0}, ray[13]={time=1,az=1},
        // ... ray[15]={time=1,az=5}
        //      k1- array of begin indx (0,12), k2- array of end indx (11,15)
        for (int i = 0; i < r.length - 1; i++) {
            time1 = r[i].getTime();
            time2 = r[i + 1].getTime();
            if (time1 != time2) {
                k2[j] = i;
                j = j + 1;
                k1[j] = i + 1;
            }
        }
        if (k2[j] < r.length - 1) {
            k1[j] = k2[j - 1] + 1;
            k2[j] = r.length - 1;
            n = j + 1;
        }

        //      if different groups have the same value of "time" (may be 2
        // and more groups) -
        //      it1= indx of "k1" of 1st group, it2= indx of "k2" of last group
        int it1 = 0, it2 = 0;
        for (int ii = 0; ii < j + 1; ii++) {
            n1 = k1[ii];
            for (int i = 0; i < j + 1; i++) {
                if (i != ii) {
                    n2 = k1[i];
                    if (r[n1].getTime() == r[n2].getTime()) {
                        it1 = ii;
                        it2 = i;
                    }
                }
            }
        }

        n1 = k1[it1];
        n2 = k1[it2];
        int s1 = k2[it1] - k1[it1] + 1;
        int s2 = k2[it2] - k1[it2] + 1;
        float[] t0 = new float[s1];
        float[] t00 = new float[s2];
        for (int i = 0; i < s1; i++) {
            t0[i] = r[n1 + i].getAz();
        }
        for (int i = 0; i < s2; i++) {
            t00[i] = r[n2 + i].getAz();
        }
        float mx0 = t0[0];
        for (int i = 0; i < s1; i++) {
            if (mx0 < t0[i]) mx0 = t0[i];
        }
        float mx00 = t00[0];
        for (int i = 0; i < s2; i++) {
            if (mx00 < t00[i]) mx00 = t00[i];
        }
        if ((mx0 > 330.0f & mx00 < 50.0f)) {
            for (int i = 0; i < s1; i++) {
                float q = r[n1 + i].getAz();
                r[n1 + i].setAz(q - 360.0f);
            }
        }
        Arrays.sort(r, new RayComparator());
        for (int i = 0; i < r.length; i++) {
            float a = r[i].getAz();
            if (a < 0 & a > -361.0f) {
                float qa = r[i].getAz();
                r[i].setAz(qa + 360.0f);
            }
        }

    }

    //  -----------------------------------------------------------------------

    /**
     * Read a date/time struct from a SIGMET file
     * @param raf file to read from
     * @return CalendarDate corresponding to the information in file
     * @throws IOException
     */
    static CalendarDate readDate(RandomAccessFile raf) throws IOException
    {
        int seconds = raf.readInt();
        int ms = raf.readShort();
        boolean is_dst = ((ms >> 10) & 0x1) == 1;
        boolean is_utc = ((ms >> 10) & 0x2) == 1;
        boolean is_local_dst = ((ms >> 10) & 0x4) == 1;
        ms = ms & 0x03FF;
        short year = raf.readShort();
        short month = raf.readShort();
        short day = raf.readShort();
        int hours = seconds / 3600;
        int minutes = (seconds / 60 - 60 * hours);
        seconds = seconds - hours * 3600 - minutes * 60;
        CalendarDate base_dt = CalendarDate.of(null, year, month, day,
                hours, minutes, seconds);
        return CalendarDate.of(base_dt.getMillis() + ms);
    }

    /**
     * Calculate azimuth of a ray
     *
     * @param az0 azimuth at beginning of ray (binary angle)
     * @param az1 azimuth at end of ray (binary angle)
     * @return float azimuth in degrees with precision of two decimal
     */
    static float calcAz(short az0, short az1) {
        // output in deg
        float azim0 = calcAngle(az0);
        float azim1 = calcAngle(az1);
        float d = Math.abs(azim0 - azim1);
        if ((az0 < 0) & (az1 > 0)) {
            d = Math.abs(360.0f - azim0) + Math.abs(azim1);
        }
        double temp = azim0 + d * 0.5;
        if (temp > 360.0) {
            temp -= 360.0;
        }
        BigDecimal bd = new BigDecimal(temp);
        BigDecimal result = bd.setScale(2, RoundingMode.HALF_DOWN);
        return result.floatValue();
    }

    /**
     * Convert 2 bytes binary angle to float
     *
     * @param angle two bytes binary angle
     * @return float value of angle in degrees with precision of two decimal
     */
    static float calcAngle(short angle) {
        final double maxval = 65536.0;
        double ang = (double) angle;
        if (ang < 0.0) {
            ang = maxval + ang;
        }
        double temp = (ang / maxval) * 360.0;
        BigDecimal bd = new BigDecimal(temp);
        BigDecimal result = bd.setScale(2, RoundingMode.HALF_DOWN);
        return result.floatValue();
    }

    /**
     * Convert 4 bytes binary angle to float
     *
     * @param ang four bytes binary angle
     * @return float value of angle with precision of two decimal in degrees
     */
    static float calcAngle(int ang) {
        final double maxval = 4294967296.0;
        double temp = (ang / maxval) * 360.0;
        BigDecimal bd = new BigDecimal(temp);
        BigDecimal result = bd.setScale(3, RoundingMode.HALF_DOWN);
        return result.floatValue();
    }

    /**
     * Calculate of Nyquist velocity
     *
     * @param prf  PRF in Hertz
     * @param wave wavelength in 1/100 of centimeters
     * @return float value of Nyquist velocity in m/sec with precision of two decimal
     */
    static float calcNyquist(int prf, int wave) {
        double tmp = (prf * wave * 0.01) * 0.25;
        tmp = tmp * 0.01;                    //Make it m/sec
        BigDecimal bd = new BigDecimal(tmp);
        BigDecimal result = bd.setScale(2, RoundingMode.HALF_DOWN);
        return result.floatValue();
    }

    class RayComparator implements Comparator<Ray> {
        public int compare(Ray ray1, Ray ray2) {
            if (ray1.getTime() < ray2.getTime()) {
                return -1;
            }
            else if (ray1.getTime() == ray2.getTime()) {
                if (ray1.getAz() < ray2.getAz()) {
                    return -1;
                }
                if (ray1.getAz() > ray2.getAz()) {
                    return 1;
                }
                if (ray1.getAz() == ray2.getAz()) {
                    return 0;
                }
            }
            else if (ray1.getTime() > ray2.getTime()) {
                return 1;
            }
            return 0;
        }
    } // class RayComparator end ----------------------------------

}