package ucar.nc2.iosp.newnids;

import ucar.nc2.Variable;

import java.util.Map;

/**
 * Created by rmay on 2/24/14.
 */
public interface Packet {
    Map<String, Variable> getData();
}
