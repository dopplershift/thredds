package ucar.nc2.iosp.newnids;

import java.util.Map;

/**
 * Created by rmay on 2/24/14.
 */
public interface ProductType {
    String name();

    String description();

    Map<String, Object> fields();
}
