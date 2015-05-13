/*
 * Copyright 1998-2014 University Corporation for Atmospheric Research/Unidata
 *
 *   Portions of this software were developed by the Unidata Program at the
 *   University Corporation for Atmospheric Research.
 *
 *   Access and use of this software shall impose the following obligations
 *   and understandings on the user. The user is granted the right, without
 *   any fee or cost, to use, copy, modify, alter, enhance and distribute
 *   this software, and any derivative works thereof, and its supporting
 *   documentation for any purpose whatsoever, provided that this entire
 *   notice appears in all copies of the software, derivative works and
 *   supporting documentation.  Further, UCAR requests that the user credit
 *   UCAR/Unidata in any publications that result from the use of this
 *   software or in any product that includes this software. The names UCAR
 *   and/or Unidata, however, may not be used in any advertising or publicity
 *   to endorse or promote any products or commercial entity unless specific
 *   written permission is obtained from UCAR/Unidata. The user also
 *   understands that UCAR/Unidata is not obligated to provide the user with
 *   any support, consulting, training or assistance of any kind with regard
 *   to the use, operation and performance of this software nor to provide
 *   the user with any updates, revisions, new versions or "bug fixes."
 *
 *   THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY EXPRESS OR
 *   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL,
 *   INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
 *   FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 *   NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
 *   WITH THE ACCESS, USE OR PERFORMANCE OF THIS SOFTWARE.
 */

package ucar.unidata.geoloc.vertical;

import ucar.ma2.*;
import ucar.ma2.ArrayDouble.D1;

import ucar.nc2.Dimension;
import ucar.nc2.Variable;
import ucar.nc2.NetcdfFile;
import ucar.unidata.util.Parameter;

import java.io.IOException;
import java.util.List;


/**
 * This implements a VerticalTransform using an existing 3D variable.
 * This is a common case when the 3D pressure or height field is stored in the file.
 *
 * @author john
 */
public class VTfromExistingData extends VerticalTransformImpl {

  /**
   * The name of the Parameter whose value is the variable that contains the 2D Height or Pressure field
   */
  public final static String existingDataField = "existingDataField";

  /**
   * The variable that contains the 2D Height or Pressure field
   */
  private Variable existingData;

  /**
   * Constructor.
   *
   * @param ds      containing Dataset
   * @param timeDim time Dimension
   * @param params  list of transformation Parameters
   */
  public VTfromExistingData(NetcdfFile ds, Dimension timeDim, List<Parameter> params) {
    super(timeDim);
    String vname = getParameterStringValue(params, existingDataField);
    this.existingData = ds.findVariable(vname);
    units = existingData.getUnitsString();
  }

  public ArrayDouble.D3 getCoordinateArray(int timeIndex) throws IOException, InvalidRangeException {
    Array data = readArray(existingData, timeIndex);

    // copy for now - better to just return Array, with promise its rank 3
    int[] shape = data.getShape();
    ArrayDouble.D3 ddata = (ArrayDouble.D3) Array.factory(DataType.DOUBLE, shape);
    MAMath.copyDouble(ddata, data);
    return ddata;
  }
  
  /**
   * Get the 1D vertical coordinate array for this time step and point
   * 
   * @param timeIndex the time index. Ignored if !isTimeDependent().
   * @param xIndex    the x index
   * @param yIndex    the y index
   * @return vertical coordinate array
   * @throws java.io.IOException problem reading data
   * @throws ucar.ma2.InvalidRangeException _more_ 
   */  
  public D1 getCoordinateArray1D(int timeIndex, int xIndex, int yIndex)
  		throws IOException, InvalidRangeException {
	  
	    
	    ArrayDouble.D3 ddata = getCoordinateArray(timeIndex);
	    int[] origin = new int[]{0, yIndex, xIndex };
	    int[] shape = new int[]{ ddata.getShape()[0] , 1, 1 };
	    
	    return (ArrayDouble.D1)ddata.section(origin, shape).reduce();	    
  }
  
}

