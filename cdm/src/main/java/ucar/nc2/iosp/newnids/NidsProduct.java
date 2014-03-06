package ucar.nc2.iosp.newnids;

import ucar.nc2.Variable;

import java.util.Map;

/**
 * Created by rmay on 2/24/14.
 */
public class NidsProduct {
    private NidsHeader header;
    private ProductDescription desc;
    private SymbologyBlock symBlock;
    private GraphicBlock graphBlock;
    private TextBlock textBlock;

    public Variable data() {
        return null;
    }

    public Map<String, Object> metadata() {
        return null;
    }
}
