package ucar.nc2.iosp.newnids;

import thredds.catalog.DataFormatType;
import ucar.ma2.Array;
import ucar.ma2.InvalidRangeException;
import ucar.ma2.Section;
import ucar.nc2.Variable;
import ucar.nc2.iosp.AbstractIOServiceProvider;
import ucar.unidata.io.RandomAccessFile;

import java.io.IOException;

/**
 * Created by rmay on 1/29/14.
 */
public class NidsIosp2 extends AbstractIOServiceProvider {
    private NidsProduct[] files;

    @Override
    public boolean isValidFile(RandomAccessFile raf) throws IOException {
        return true;
    }

    @Override
    public Array readData(Variable v2, Section section) throws IOException,
            InvalidRangeException {
        return null;
    }

    @Override
    public String getFileTypeId() {
        return DataFormatType.NIDS.toString();
    }

    @Override
    public String getFileTypeDescription() {
        return "NEXRAD Level-III (NIDS) Products";
    }
}
