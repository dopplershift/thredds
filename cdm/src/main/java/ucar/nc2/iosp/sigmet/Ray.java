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

import ucar.ma2.IndexIterator;
import ucar.ma2.Range;
import ucar.unidata.io.RandomAccessFile;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Comparator;
import java.util.Formatter;
import java.util.Map;

/**
 * @author yuanho
 * @since Apr 7, 2010
 */
public class Ray {
    private short bins;
    int dataRead;
    int offset;
    int offset1;
    private float range, step, az, elev;
    private short time;
    //  private float[] val;
    int nsweep;
    int datatype;

    public Ray(short bins, int nsweep, int datatype) {
        this(-999.99f, -999.99f, -999.99f, -999.99f, bins, (short) (-99),
                -999, 0, -999, nsweep, datatype);
    }

    public Ray(float range, float step, float az, float elev, short bins,
               short time, int offset, int dataRead,
               int offset1, int nsweep, int datatype) {

        //  this.val = new float[bins];
        setRange(range);
        setStep(step);
        setAz(az);
        setElev(elev);
        setBins(bins);
        setTime(time);
        setOffset(offset);
        setDataRead(dataRead);
        setOffset1(offset1);
        // setVal(val);
        setNsweep(nsweep);
        setDataType(datatype);
    }

    public int getDataType() {
        return datatype;
    }

    public void setDataType(int datatype) {
        this.datatype = datatype;
    }

    public float getRange() {
        return range;
    }

    public void setRange(float range) {
        this.range = range;
    }

    public float getStep() {
        return step;
    }

    public void setStep(float step) {
        this.step = step;
    }

    public int getNsweep() {
        return nsweep;
    }

    public void setNsweep(int nsweep) {
        this.nsweep = nsweep;
    }

    public float getAz() {
        if (az < 0.0f & az > -361.0f) {
            az = 360.0f + az;
        }

        return az;
    }

    public void setAz(float az) {
        this.az = az;
    }

    public float getElev() {
        return elev;
    }

    public void setElev(float elev) {
        this.elev = elev;
    }

    public short getBins() {
        return bins;
    }

    public void setBins(short bins) {
        this.bins = bins;
    }

    public short getTime() {
        return time;
    }

    public void setTime(short time) {
        this.time = time;
    }

    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }

    public int getDataRead() {
        return dataRead;
    }

    public void setDataRead(int dataRead) {
        this.dataRead = dataRead;
    }

    public int getOffset1() {
        return offset1;
    }

    public void setOffset1(int offset1) {
        this.offset1 = offset1;
    }

    /*   public float[] getVal() {
          return val;
      }

      public void setVal(float[] val) {
          System.arraycopy(val, 0, this.val, 0, bins);
      }
    */

    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        else if (o instanceof Ray) {
            Ray oo = (Ray) o;

            return (range == oo.range & step == oo.step & az == oo.az & elev
                    == oo.elev & bins == oo.bins
                    & time == oo.time);
        }
        else {
            return false;
        }
    }

    public int hashCode() {
        return new Float(range).hashCode() + new Float(step).hashCode() + new
                Float(az).hashCode()
                + new Float(elev).hashCode() + new Short(bins).hashCode() +
                new Short(time).hashCode();

        // val.hashCode();
    }

    public String toString() {
        Formatter sb = new Formatter();

        sb.format("Range=%f Step=%f", range, step);
        if (az > -361.0f & az < 0.0f) {
            az = 360.0f + az;
        }

        sb.format(" Az=%f Elev=%f Bins=%d Time=%d", az, elev, bins, time);
        // for (int i=0; i<bins; i++) { sb.append(" "+val[i]); }
        return sb.toString();
    }


    /**
     * Read data from this record.
     *
     * @param raf       read from this file
     * @param gateRange handles the possible subset of data to return
     * @param ii        put the data here
     * @throws java.io.IOException on read error
     */
    public void readData(RandomAccessFile raf,
                         SigmetVolumeScan.RecordsHeader recHdr, Range gateRange,
                         IndexIterator ii) throws IOException {
        final int REC_SIZE = 6144;
        raf.seek(offset);
        byte[] data = new byte[bins];
        float[] dd = new float[bins];
        byte d;
        int nb = 0;
        short a00;
        // raf.readFully(data);
        if (dataRead > 0) {
            raf.seek(offset);

            for (int i = 0; i < dataRead; i++) {
                d = raf.readByte();
                dd[i] = calcData(getDataType(), d, recHdr.velScale);
                nb++;
            }
        }
        // System.out.println("this is az " + getAz());
        raf.seek(offset1);
        int cur_len = offset1;

        while (nb < (int) bins) {
            // --- Check if the code=1 ("1" means an end of a ray)
            a00 = raf.readShort();
            cur_len = cur_len + 2;

            if (a00 == (short) 1) {
                for (int uk = 0; uk < (int) bins; uk++) {
                    dd[uk] = -999.99f;
                }
                break;
            }

            if (a00 < 0) {    // -- This is data
                int nwords = a00 & 0x7fff;
                int dataRead1 = nwords * 2;
                int pos = 0;
                if (cur_len % REC_SIZE == 0) {
                    pos = 0;
                    break;
                }
                raf.seek(cur_len);
                for (int i = 0; i < dataRead1; i++) {
                    d = raf.readByte();
                    dd[nb] = calcData(getDataType(), d, recHdr.velScale);
                    nb = nb + 1;
                    cur_len = cur_len + 1;
                    if (nb % REC_SIZE == 0) {
                        pos = i + 1;
                        break;
                    }
                }

                if (pos > 0) {
                    break;
                }
            }
            else if (a00 > 0 & a00 != 1) {
                int num_zero = a00 * 2;
                int dataRead1 = num_zero;
                for (int k = 0; k < dataRead1; k++) {
                    dd[nb + k] = calcData(getDataType(), (byte) 0, recHdr.velScale);
                }
                nb = nb + dataRead1;
                if (cur_len % REC_SIZE == 0) {
                    break;
                }

            }
        }                     // ------ end of while for
        // num_bins---------------------------------

        for (int i = gateRange.first(); i <= gateRange.last(); i += gateRange
                .stride()) {
            if (i >= bins)
                ii.setFloatNext(Float.NaN);
            else
                ii.setFloatNext(dd[i]);
        }

    } // end of readData

    /**
     * Calculate data values from raw ingest data
     *
     * @param dty    type of data ( "Total_Power", "Reflectivity", "Velocity",
     *               "Width", "Differential_Reflectivity")
     * @param data   1-byte input value
     * @param velScale scaling value for velocity data
     * @return float value with precision of two decimal
     */
    static float calcData(int dty, byte data, float velScale) {
        double temp = -999.99;
        switch (dty) {
            default:        // dty=1,2 -total_power, reflectivity (dBZ)
                if (data != 0) {
                    temp = (((int) data & 0xFF) - 64) * 0.5;
                }
                break;
            case 3:        // dty=3 - mean velocity (m/sec)
                if (data != 0) {
                    temp = ((((int) data & 0xFF) - 128) / 127.0) * velScale;
                }
                break;
            case 4:        // dty=4 - spectrum width (m/sec)
                if (data != 0) {
                    double v = ((((int) data & 0xFF) - 128) / 127.0) * velScale;
                    temp = (((int) data & 0xFF) / 256.0) * v;
                }
                break;
            case 5:        // dty=5 - differential reflectivity (dB)
                if (data != 0) {
                    temp = ((((int) data & 0xFF) - 128) / 16.0);
                }
                break;
        }
        BigDecimal bd = new BigDecimal(temp);
        BigDecimal result = bd.setScale(2, RoundingMode.HALF_DOWN);
        return result.floatValue();
    }

}    // class Ray end------------------------------------------
