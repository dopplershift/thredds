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
    public List<List<Float>> azimuth, elevation;
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

                raf.skipBytes(2);
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
        azimuth = new ArrayList<>();
        elevation = new ArrayList<>();
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
            while (azimuth.size() < nsweep)
            {
                azimuth.add(new ArrayList<Float>());
                elevation.add(new ArrayList<Float>());
            }
            azimuth.get(nsweep - 1).add(az);
            elevation.get(nsweep - 1).add(elev);
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
                    ray = new Ray(num_bins, nsweep, dty);
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

//        groups.put("Times", sortScans(time, 1000));
    }

    private void addRay(int dataType, Ray ray)
    {
        // Convert data type Id to name
        String datStr = dataTypeNames[dataType];

        // Get list of sweeps for this type or create if necessary
        List<List<Ray>> data = groups.get(datStr);
        if (data == null)
        {
            data = new ArrayList<>(this.recHdr.num_sweeps);
        }

        while (data.size() <= ray.getNsweep())
            data.add(new ArrayList<Ray>(this.recHdr.num_rays));

        List<Ray> sweep = data.get(ray.getNsweep());
        // TODO: Need to handle missing rays by actually getting ray index
        // from file.
//        while (sweep.size() <= ray.getOffset())
//            sweep.add(new Ray(ray.getBins(), ray.getNsweep(),
//                    ray.getDataType()));

        sweep.add(ray);
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
}