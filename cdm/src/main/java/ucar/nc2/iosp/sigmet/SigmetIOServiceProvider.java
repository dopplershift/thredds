package ucar.nc2.iosp.sigmet;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.ByteBuffer;
import java.nio.channels.WritableByteChannel;
import java.util.ArrayList;
import java.util.List;

import ucar.ma2.Range;
import ucar.ma2.Array;
import ucar.ma2.ArrayFloat;
import ucar.ma2.ArrayInt;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.ma2.IndexIterator;
import ucar.ma2.InvalidRangeException;
import ucar.ma2.Section;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;
import ucar.nc2.constants.AxisType;
import ucar.nc2.constants.CDM;
import ucar.nc2.constants._Coordinate;
import ucar.nc2.iosp.AbstractIOServiceProvider;
import ucar.nc2.iosp.Layout;
import ucar.nc2.iosp.LayoutRegular;
import ucar.nc2.time.CalendarDate;
import ucar.nc2.time.CalendarDateUnit;
import ucar.nc2.time.CalendarPeriod;
import ucar.unidata.io.RandomAccessFile;

/**
 * Implementation of the ServerProvider pattern provides input/output
 * of the SIGMET-IRIS dataset. IOSP are managed by the NetcdfFile class.
 * When the SigmetDataset is requested by calling NetcdfFile.open(), the file
 * is opened as a ucar.unidata.io.RandomAccessFile.
 * The SIGMET-IRIS data format are described in "IRIS Programmer's Manual" ch.3
 * The SIGMET-IRIS file consists of records with fixed length=6144 bytes.
 * Data is written in
 * little-endian format. The organization of raw product SIGMET-IRIS file is:
 * Record #1   { <product_hdr> 0,0,0,...}
 * Record #2   { <ingest_header> 0,0,0,... }
 * ---if there are several sweeps (usually 24) and one type of data:
 * Record #3   { <raw_prod_bhdr><ingest_data_header>...Data for sweep 1.. }
 * Record #4   { <raw_prod_bhdr>...Data for sweep 1...  }
 * .............................................
 * Record #n   { <raw_prod_bhdr>...Data for sweep 1... 0.... }
 * Record #n+1 { <raw_prod_bhdr><ingest_data_header>...Data for sweep 2.. }
 * Record #n+2 { <raw_prod_bhdr>...Data for sweep 2...  }
 * .............................................
 * Record #m   { <raw_prod_bhdr>...Data for sweep 2... 0...  }
 * Record #m+1 { <raw_prod_bhdr><ingest_data_header>...Data for sweep 3.. }
 * .............................................
 * Structure of "Data for sweep" is:  <ray_header><ray_data>..
 * .<ray_header><ray_data>...
 * <ray_header> and <ray_data> are encoded with the compression algorithm
 * ("IRIS Programmer's Manual" 3.5.4.1)
 * ---if there are "n" types of data (usually 4 or 5) and one sweep:
 * Record #3   { <raw_prod_bhdr><ingest_data_header(data_type_1)
 * ><ingest_data_header(data_type_2)>...
 * <ingest_data_header(data_type_n)>...Data...}
 * Record #4   { <raw_prod_bhdr>...Data...  }
 * .............................................
 * Record #n   { <raw_prod_bhdr>...Data...  }
 * Structure of "Data" is:  <ray_header(data_type_1)><ray_data(data_type_1)
 * ><ray_header(data_type_2)><ray_data(data_type_2)>...
 * <ray_header(data_type_n)><ray_data(data_type_n)><ray_header(data_type_1)
 * ><ray_data(data_type_1)>
 * <ray_header(data_type_2)><ray_data(data_type_2)>... <ray_header
 * (data_type_n)><ray_data(data_type_n)>...
 * <ray_header> and <ray_data> are encoded with the compression algorithm
 * ("IRIS Programmer's Manual" 3.5.4.1)
 *
 * @author yuanho
 * @see "ftp://ftp.sigmet.com/outgoing/manuals/IRIS_Programmers_Manual.pdf"
 * esp Chapter 4
 */

public class SigmetIOServiceProvider extends AbstractIOServiceProvider {
    private ArrayList<Variable> varList = null;
    private int[] sweep_bins = null;

    private SigmetVolumeScan volScan;

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

    public static void main(String[] args) {
        String infile = " ";
        if (args.length == 1) {
            infile = args[0];
        }
        else {
            System.out.println("Usage: java SigmetIOServiceProvider inputFile");
            System.exit(0);
        }
        try {
            NetcdfFile.registerIOProvider(SigmetIOServiceProvider.class);
            NetcdfFile ncfile = NetcdfFile.open(infile);
            System.out.println("ncfile = \n" + ncfile);
        } catch (Exception e) {
            System.out.println("MAIN!!!   " + e.toString());
            e.printStackTrace();
        }
    }

    public String getFileTypeDescription() {
        return "SIGMET-IRIS";
    }

    public String getFileTypeVersion() {
        return "SIGMET-IRIS";
    }

    public String getFileTypeId() {
        return "SIGMET";
    }

    /**
     * Check if this is a valid SIGMET-IRIS file for this IOServiceProvider.
     */
    public boolean isValidFile(ucar.unidata.io.RandomAccessFile raf) {
        try {
            raf.order(RandomAccessFile.LITTLE_ENDIAN);
            // The first struct in the file is the product_hdr,
            // which will have the
            // standard structure_header, followed by other embedded structures.
            // Each of these structures also have a structure header. To
            // validate
            // the file we check for a product_hdr (by looking for type 27 in
            // the
            // structure_header), then a product_configuration structure (by
            // looking
            // for type 26 in its structure_header), then checking that that
            // the product_configuration does indicate a type of RAW data
            // (type 15)
            raf.seek(0);
            short[] data = new short[13];
            raf.readShort(data, 0, 13);
            return (data[0] == (short) 27 &&
                    data[6] == (short) 26 &&
                    data[12] == (short) 15);
        } catch (IOException ioe) {
            System.out.println("In isValidFile(): " + ioe.toString());
            return false;
        }
    }

    /**
     * Open existing file, and populate ncfile with it.
     */
    public void open(ucar.unidata.io.RandomAccessFile raf,
                     ucar.nc2.NetcdfFile ncfile,
                     ucar.nc2.util.CancelTask cancelTask) throws java.io
            .IOException {
        super.open(raf, ncfile, cancelTask);
        volScan = new SigmetVolumeScan(raf, ncfile);
        varList = init(raf, ncfile);
    }

    /**
     * Define Dimensions, Variables, Attributes in ncfile
     *
     * @param raf    ucar.unidata.io.RandomAccessFile corresponds of SIGMET
     *               datafile.
     * @param ncfile an empty NetcdfFile object which will be filled.
     * @return ArrayList of Variables of ncfile
     */
    public ArrayList<Variable> init(ucar.unidata.io.RandomAccessFile raf,
                                    ucar.nc2.NetcdfFile ncfile) throws java
            .io.IOException {
        SigmetVolumeScan.RecordsHeader recHdr = volScan.header();

        short num_rays = recHdr.num_rays;
        short bins = recHdr.num_bins;
        float range_first = recHdr.range_first;
        float range_last = recHdr.range_last;
        short number_sweeps = recHdr.num_sweeps;

        // define number of gates for every sweep
        sweep_bins = volScan.getNumberGates();

        // add simple dimensions for scan (sweep) and radial
        Dimension scanR = new Dimension("sweep", number_sweeps, true);
        ncfile.addDimension(null, scanR);

        Dimension radial = new Dimension("radial", num_rays, true);
        ncfile.addDimension(null, radial);

        // Add a gates dimension that is sized to the maximum number
        int max_bins = sweep_bins[0];
        for (int j = 1; j < number_sweeps; j++) {
            max_bins = java.lang.Math.max(max_bins, sweep_bins[j]);
        }
        Dimension gateR = new Dimension("gates", max_bins, true);
        ncfile.addDimension(null, gateR);

        // Set up the dimensions for the data fields
        ArrayList<Dimension> dims = new ArrayList<>();
        dims.add(scanR);
        dims.add(radial);
        dims.add(gateR);

        // Storage for variables.
        ArrayList<Variable> varList = new ArrayList<>();

        // Create a variable for each parameter in the file
        int nparams = recHdr.params.size();
        String coordinates = "time elevationR azimuthR distanceR";
        for (int param: recHdr.params) {
            ++nparams;
            String var_name = dataTypeNames[param];
            Variable var = new Variable(ncfile, null, null, var_name);
            var.setDataType(DataType.FLOAT);
            var.setDimensions(dims);
            var.addAttribute(new Attribute(CDM.LONG_NAME, var_name));
            var.addAttribute(new Attribute(CDM.UNITS,
                    dataTypeUnit[param]));
            var.addAttribute(new Attribute(_Coordinate.Axes, coordinates));
            var.addAttribute(new Attribute(CDM.MISSING_VALUE, -999.99f));
            ncfile.addVariable(null, var);
            varList.add(var);
        }

        CalendarDateUnit time_unit = CalendarDateUnit.of(
                volScan.start_date.getCalendar(), CalendarPeriod.Field.Millisec,
                volScan.start_date);

        dims.clear();
        dims.add(scanR);
        dims.add(radial);

        // add "time" variable
        Variable time = new Variable(ncfile, null, null, "time");
        time.setDataType(DataType.INT);
        time.setDimensions(dims);
        time.addAttribute(new Attribute(CDM.LONG_NAME,
                "time from start of volume"));
        time.addAttribute(new Attribute(CDM.UNITS, time_unit.toString()));
        time.addAttribute(new Attribute(_Coordinate.AxisType,
                AxisType.Time.toString()));
        time.addAttribute(new Attribute(CDM.MISSING_VALUE, -99));
        ncfile.addVariable(null, time);
        varList.add(time);

        // add "elevationR" variable
        Variable elevationR = new Variable(ncfile, null, null, "elevationR");
        elevationR.setDataType(DataType.FLOAT);
        elevationR.setDimensions(dims);
        elevationR.addAttribute(new Attribute(CDM.LONG_NAME,
                "elevation angle"));
        elevationR.addAttribute(new Attribute(CDM.UNITS, "degrees"));
        elevationR.addAttribute(new Attribute(_Coordinate.AxisType,
                AxisType.RadialElevation.toString()));
        elevationR.addAttribute(new Attribute(CDM.MISSING_VALUE, -999.99f));
        ncfile.addVariable(null, elevationR);
        varList.add(elevationR);

        // add "azimuthR" variable
        Variable azimuthR = new Variable(ncfile, null, null, "azimuthR");
        azimuthR.setDataType(DataType.FLOAT);
        azimuthR.setDimensions(dims);
        azimuthR.addAttribute(new Attribute(CDM.LONG_NAME, "azimuth angle"));
        azimuthR.addAttribute(new Attribute(CDM.UNITS, "degrees"));
        azimuthR.addAttribute(new Attribute(_Coordinate.AxisType,
                AxisType.RadialAzimuth.toString()));
        azimuthR.addAttribute(new Attribute(CDM.MISSING_VALUE, -999.99f));
        ncfile.addVariable(null, azimuthR);
        varList.add(azimuthR);

        // add "distanceR" variable. This is a function of sweep and gate
        dims.clear();
        dims.add(scanR);
        dims.add(gateR);
        Variable distanceR = new Variable(ncfile, null, null, "distanceR");
        distanceR.setDataType(DataType.FLOAT);
        distanceR.setDimensions(dims);
        distanceR.addAttribute(new Attribute(CDM.LONG_NAME, "radial distance"));
        distanceR.addAttribute(new Attribute(CDM.UNITS, "m"));
        distanceR.addAttribute(new Attribute(_Coordinate.AxisType,
                AxisType.RadialDistance.toString()));
        ncfile.addVariable(null, distanceR);
        varList.add(distanceR);

        // add "numGates" variable
        dims.clear();
        dims.add(scanR);
        Variable numGates = new Variable(ncfile, null, null, "numGates");
        numGates.setDataType(DataType.INT);
        numGates.setDimensions(dims);
        numGates.addAttribute(new Attribute(CDM.LONG_NAME,
                "number of gates in the sweep"));
        ncfile.addVariable(null, numGates);
        varList.add(numGates);

        // add global attributes
        ncfile.addAttribute(null, new Attribute("definition",
                "SIGMET-IRIS RAW"));
        ncfile.addAttribute(null, new Attribute("description",
                "SIGMET-IRIS data are read by NetCDF IOSP"));
        ncfile.addAttribute(null, new Attribute("StationName",
                recHdr.stnName));
        ncfile.addAttribute(null, new Attribute("StationName_SetupUtility",
                recHdr.stnName_util));
        ncfile.addAttribute(null, new Attribute("radar_lat",
                recHdr.radar_lat));
        ncfile.addAttribute(null, new Attribute("radar_lon",
                recHdr.radar_lon));
        ncfile.addAttribute(null, new Attribute("ground_height",
                recHdr.ground_height));
        ncfile.addAttribute(null, new Attribute("radar_height",
                recHdr.radar_height));
        ncfile.addAttribute(null, new Attribute("radar_alt",
                recHdr.radar_alt / 100));
        ncfile.addAttribute(null, new Attribute("time_coverage_start",
                volScan.start_date.toString()));
        ncfile.addAttribute(null, new Attribute("time_coverage_end",
                volScan.end_date.toString()));
        ncfile.addAttribute(null, new Attribute("num_data_types", nparams));
        ncfile.addAttribute(null, new Attribute("number_sweeps",
                number_sweeps));
        ncfile.addAttribute(null, new Attribute("num_rays", num_rays));
        ncfile.addAttribute(null, new Attribute("max_number_gates", bins));
        ncfile.addAttribute(null, new Attribute("range_first", range_first));
        ncfile.addAttribute(null, new Attribute("range_last", range_last));
        ncfile.addAttribute(null, new Attribute("DataType", "Radial"));
        ncfile.addAttribute(null, new Attribute(CDM.CONVENTIONS,
                _Coordinate.Convention));

        // --------- fill all of values in the ncfile ------
        doNetcdfFileCoordinate(ncfile, varList, recHdr);

        ncfile.finish();

        return varList;
    }

    /**
     * Fill all of the variables/attributes in the ncfile
     *
     * @param ncfile  NetcdfFile object which will be filled.
     * @param varList ArrayList of Variables of ncfile
     * @param recHdr  java.util.Map with values for Attributes
     */
    private void doNetcdfFileCoordinate(ucar.nc2.NetcdfFile ncfile,
                                       ArrayList<Variable> varList,
                                       SigmetVolumeScan.RecordsHeader recHdr) {
        // prepare attribute values

        short num_rays = recHdr.num_rays;
        float range_first = recHdr.range_first;
        float range_last = recHdr.range_last;
        short number_sweeps = recHdr.num_sweeps;

        // set all of Variables
        try {
            Ray[] rtemp = new Ray[(int) num_rays];

            Variable[] distanceR = new Variable[number_sweeps];
            ArrayFloat.D1[] distArr = new ArrayFloat.D1[number_sweeps];
            Index[] distIndex = new Index[number_sweeps];
            String distName = "distanceR";
            for (int i = 0; i < number_sweeps; i++) {
                if (number_sweeps > 1) {
                    distName = "distanceR_sweep_" + (i + 1);
                }
                for (Variable aVarList : varList) {
                    if ((aVarList.getShortName()).equals(distName.trim())) {
                        distanceR[i] = aVarList;
                        break;
                    }
                }
                distArr[i] = (ArrayFloat.D1) Array.factory(DataType.FLOAT,
                        distanceR[i].getShape());
                distIndex[i] = distArr[i].getIndex();

                int ngates = sweep_bins[i];
                float stp = calcStep(range_first, range_last, (short) ngates);
                for (int ii = 0; ii < ngates; ii++) {
                    distArr[i].setFloat(distIndex[i].set(ii),
                            (range_first + ii * stp));
                }
            }
            List rgp = volScan.getTotalPowerGroups();
            if (rgp.size() == 0) rgp = volScan.getReflectivityGroups();
            List[] sgp = new ArrayList[number_sweeps];
            for (int i = 0; i < number_sweeps; i++) {
                sgp[i] = (List) rgp.get((short) i);
            }


            Variable[] time = new Variable[number_sweeps];
            ArrayInt.D1[] timeArr = new ArrayInt.D1[number_sweeps];
            Index[] timeIndex = new Index[number_sweeps];
            String t_n = "time";
            for (int i = 0; i < number_sweeps; i++) {
                if (number_sweeps > 1) {
                    t_n = "time_sweep_" + (i + 1);
                }
                for (Variable aVarList : varList) {
                    if ((aVarList.getShortName()).equals(t_n.trim())) {
                        time[i] = aVarList;
                        break;
                    }
                }

                //                if (time[i].getShape().length == 0) {
                //                    continue;
                //                }
                timeArr[i] = (ArrayInt.D1) Array.factory(DataType.INT,
                        time[i].getShape());
                timeIndex[i] = timeArr[i].getIndex();
                List rlist = sgp[i];

                for (int jj = 0; jj < num_rays; jj++) {
                    rtemp[jj] = (Ray) rlist.get(jj);
                }    //ray[i][jj]; }
                for (int jj = 0; jj < num_rays; jj++) {
                    timeArr[i].setInt(timeIndex[i].set(jj),
                            rtemp[jj].getTime());
                }
            }

            // NCdump.printArray(timeArr[0], "time", System.out, null);

            Variable[] azimuthR = new Variable[number_sweeps];
            ArrayFloat.D1[] azimArr = new ArrayFloat.D1[number_sweeps];
            Index[] azimIndex = new Index[number_sweeps];
            String azimName = "azimuthR";
            for (int i = 0; i < number_sweeps; i++) {
                if (number_sweeps > 1) {
                    azimName = "azimuthR_sweep_" + (i + 1);
                }
                for (Variable aVarList : varList) {
                    if ((aVarList.getShortName()).equals(azimName.trim())) {
                        azimuthR[i] = aVarList;
                        break;
                    }
                }
                azimArr[i] = (ArrayFloat.D1) Array.factory(DataType.FLOAT,
                        azimuthR[i].getShape());
                azimIndex[i] = azimArr[i].getIndex();
                List rlist = sgp[i];

                for (int jj = 0; jj < num_rays; jj++) {
                    rtemp[jj] = (Ray) rlist.get(jj);
                } //ray[i][jj]; }
                for (int jj = 0; jj < num_rays; jj++) {
                    azimArr[i].setFloat(azimIndex[i].set(jj),
                            rtemp[jj].getAz());
                }
            }
            //NCdump.printArray(azimArr[0], "azimuthR", System.out, null);

            Variable[] elevationR = new Variable[number_sweeps];
            ArrayFloat.D1[] elevArr = new ArrayFloat.D1[number_sweeps];
            Index[] elevIndex = new Index[number_sweeps];
            String elevName = "elevationR";
            for (int i = 0; i < number_sweeps; i++) {
                if (number_sweeps > 1) {
                    elevName = "elevationR_sweep_" + (i + 1);
                }
                for (Variable aVarList : varList) {
                    if ((aVarList.getShortName()).equals(elevName.trim())) {
                        elevationR[i] = aVarList;
                        break;
                    }
                }
                elevArr[i] = (ArrayFloat.D1) Array.factory(DataType.FLOAT,
                        elevationR[i].getShape());
                elevIndex[i] = elevArr[i].getIndex();
                List rlist = sgp[i];

                for (int jj = 0; jj < num_rays; jj++) {
                    rtemp[jj] = (Ray) rlist.get(jj);
                } //ray[i][jj]; }
                for (int jj = 0; jj < num_rays; jj++) {
                    elevArr[i].setFloat(elevIndex[i].set(jj),
                            rtemp[jj].getElev());
                }
            }
            // NCdump.printArray(elevArr[0], "elevationR", System.out, null);

            Variable numGates = null;
            for (int i = 0; i < number_sweeps; i++) {
                for (Variable aVarList : varList) {
                    if ((aVarList.getShortName()).equals("numGates")) {
                        numGates = aVarList;
                        break;
                    }
                }
            }
            ArrayInt.D1 gatesArr = (ArrayInt.D1) Array.factory(DataType.INT,
                    numGates.getShape());
            Index gatesIndex = gatesArr.getIndex();

            for (int i = 0; i < number_sweeps; i++) {
                List rlist = sgp[i];
                for (int jj = 0; jj < num_rays; jj++) {
                    rtemp[jj] = (Ray) rlist.get(jj);
                } //ray[i][jj]; }
                gatesArr.setInt(gatesIndex.set(i), rtemp[0].getBins());
            }

            for (int i = 0; i < number_sweeps; i++) {
                distanceR[i].setCachedData(distArr[i], false);
                time[i].setCachedData(timeArr[i], false);
                azimuthR[i].setCachedData(azimArr[i], false);
                elevationR[i].setCachedData(elevArr[i], false);
            }
            numGates.setCachedData(gatesArr, false);
            // startSweep.setCachedData(sweepArr, false);

            //          -------------------------------------------------
            // int b=(int)ray[0][0].getBins();

            // -- Test of readData() and readToByteChannel() -----------------
      /*
Range r1=new Range(356, 359);
Range r2=new Range(0, 15);
java.util.List arlist=new ArrayList();
arlist.add(r1);
arlist.add(r2);
Array testArr=readData(v[0], new Section(arlist));
NCdump.printArray(testArr, "Total_Power_sweep_1", System.out, null);
WritableByteChannel channel=new FileOutputStream(new File("C:\\netcdf\\tt
.dat")).getChannel();
long ikk=readToByteChannel(v[0], new Section(arlist), channel);
System.out.println("IKK="+ikk);
channel.close();
      */
            //---------------------------------------------------

        } catch (Exception e) {
            System.out.println(e.toString());
            e.printStackTrace();
        }
    }   //----------- end of doNetcdf ----------------------------------

    /**
     * Read data from a top level Variable and return a memory resident Array.
     *
     * @param v2      Variable. It may have FLOAT/INTEGER data type.
     * @param section wanted section of data of Variable. The section list is a
     *                list
     *                of ucar.ma2.Range which define the requested data subset.
     * @return Array of data which will be read from Variable through this call.
     */
    public Array readData1(ucar.nc2.Variable v2, Section section)
            throws IOException, InvalidRangeException {
        //doData(raf, ncfile, varList);
        int[] sh = section.getShape();
        Array temp = Array.factory(v2.getDataType(), sh);
        long pos0 = 0;
        // Suppose that the data has LayoutRegular
        LayoutRegular index = new LayoutRegular(pos0, v2.getElementSize(),
                v2.getShape(), section);
        if (v2.getShortName().startsWith("time") | v2.getShortName()
                .startsWith("numGates")) {
            temp = readIntData(index, v2);
        }
        else {
            temp = readFloatData(index, v2);
        }
        return temp;
    }

    public Array readData(Variable v2, Section section) throws IOException,
            InvalidRangeException {
        // Vgroup vgroup = (Vgroup) v2.getSPobject();
        // Range scanRange = section.getRange(0);
        // Range radialRange = section.getRange(1);
        // Range gateRange = section.getRange(2);

        Array data = Array.factory(v2.getDataType().getPrimitiveClassType(),
                section.getShape());
        IndexIterator ii = data.getIndexIterator();

        List<List<Ray>> groups;
        String shortName = v2.getShortName();
        if (shortName.startsWith("Reflectivity"))
            groups = volScan.getReflectivityGroups();
        else if (shortName.startsWith("Velocity"))
            groups = volScan.getVelocityGroups();
        else if (shortName.startsWith("TotalPower"))
            groups = volScan.getTotalPowerGroups();
        else if (shortName.startsWith("Width"))
            groups = volScan.getWidthGroups();
        else if (shortName.startsWith("DiffReflectivity"))
            groups = volScan.getDifferentialReflectivityGroups();
        else
            throw new IllegalStateException("Illegal variable name = " +
                    shortName);

        if (section.getRank() == 2) {
            Range radialRange = section.getRange(0);
            Range gateRange = section.getRange(1);
            List<Ray> lli = groups.get(0);
            readOneScan(lli, radialRange, gateRange, ii);
        }
        else {
            Range scanRange = section.getRange(0);
            Range radialRange = section.getRange(1);
            Range gateRange = section.getRange(2);
            for (int i = scanRange.first(); i <= scanRange.last(); i +=
                    scanRange.stride()) {
                readOneScan(groups.get(i), radialRange, gateRange, ii);
            }
        }
        return data;
    }


    private void readOneScan(List<Ray> mapScan, Range radialRange,
                             Range gateRange, IndexIterator ii) throws
            IOException {
        int siz = mapScan.size();
        for (int i = radialRange.first(); i <= radialRange.last(); i +=
                radialRange.stride()) {
            if (i >= siz)
                readOneRadial(null, gateRange, ii);
            else {
                Ray r = mapScan.get(i);
                readOneRadial(r, gateRange, ii);
            }
        }
    }

    private void readOneRadial(Ray r, Range gateRange,
                               IndexIterator ii) throws IOException {
        if (r == null) {
            for (int i = gateRange.first(); i <= gateRange.last(); i +=
                    gateRange.stride())
                ii.setFloatNext(Float.NaN);
            return;
        }
        r.readData(volScan.raf, volScan.header(), gateRange, ii);
    }

    /**
     * Read data from a top level Variable of INTEGER data type and return a
     * memory resident Array.
     *
     * @param index LayoutRegular object
     * @param v2    Variable has INTEGER data type.
     * @return Array of data which will be read from Variable through this call.
     */
    public Array readIntData(LayoutRegular index, Variable v2)
            throws IOException {
        int[] var = (int[]) (v2.read().get1DJavaArray(v2.getDataType()
                .getPrimitiveClassType()));
        int[] data = new int[(int) index.getTotalNelems()];
        while (index.hasNext()) {
            Layout.Chunk chunk = index.next();
            System.arraycopy(var, (int) chunk.getSrcPos() / 4, data,
                    (int) chunk.getDestElem(), chunk.getNelems());
        }
        return Array.factory(data);
    }

    /**
     * Read data from a top level Variable of FLOAT data type and return a
     * memory resident Array.
     *
     * @param index LayoutRegular object
     * @param v2    Variable has FLOAT data type.
     * @return Array of data which will be read from Variable through this call.
     */
    public Array readFloatData(LayoutRegular index, Variable v2)
            throws IOException {
        float[] var = (float[]) (v2.read().get1DJavaArray(v2.getDataType()
                .getPrimitiveClassType()));
        float[] data = new float[(int) index.getTotalNelems()];
        while (index.hasNext()) {
            Layout.Chunk chunk = index.next();
            System.arraycopy(var, (int) chunk.getSrcPos() / 4, data,
                    (int) chunk.getDestElem(), chunk.getNelems());
        }
        return Array.factory(data);
    }
    //----------------------------------------------------------------------------------

    /**
     * Read data from a top level Variable and send data to a
     * WritableByteChannel.
     *
     * @param v2      Variable
     * @param section wanted section of data of Variable. The section list is a
     *                list
     *                of ucar.ma2.Range which define the requested data subset.
     * @param channel WritableByteChannel object - channel that can write bytes.
     * @return the number of bytes written, possibly zero.
     */
    public long readToByteChannel11(ucar.nc2.Variable v2, Section section, WritableByteChannel channel)
            throws java.io.IOException, ucar.ma2.InvalidRangeException {
        Array data = readData(v2, section);
        float[] ftdata = new float[(int) data.getSize()];
        byte[] bytedata = new byte[(int) data.getSize() * 4];
        IndexIterator iter = data.getIndexIterator();
        int i = 0;
        ByteBuffer buffer = ByteBuffer.allocateDirect(bytedata.length);
        while (iter.hasNext()) {
            ftdata[i] = iter.getFloatNext();
            bytedata[i] = new Float(ftdata[i]).byteValue();
            buffer.put(bytedata[i]);
            i++;
        }
        buffer = ByteBuffer.wrap(bytedata);
        // write the bytes to the channel
        int count = channel.write(buffer);
        System.out.println("COUNT=" + count);
        // check if all bytes where written
        if (buffer.hasRemaining()) {
            // if not all bytes were written, move the unwritten bytes to the beginning and
            // set position just after the last unwritten byte
            buffer.compact();
        }
        else {
            buffer.clear();
        }
        return (long) count;
    }


    /**
     * Calculate distance between sequential bins in a ray
     *
     * @param range_first range of first bin in centimeters
     * @param range_last  range of last bin in centimeters
     * @param num_bins    number of bins
     * @return float distance in centimeters with precision of two decimal
     */
    static float calcStep(float range_first, float range_last, short num_bins) {
        float step = (range_last - range_first) / (num_bins - 1);
        BigDecimal bd = new BigDecimal(step);
        BigDecimal result = bd.setScale(2, RoundingMode.HALF_DOWN);
        return result.floatValue();
    }

    /**
     * Calculate time as hh:mm:ss
     *
     * @param t  number of seconds since midnight for start of sweep
     * @param t0 time in seconds from start of sweep
     * @return time as string "hh:mm:ss"
     */
    static String calcTime(int t, int t0) {
        StringBuilder tim = new StringBuilder();
        int[] tt = new int[3];
        int mmh = (t + t0) / 60;
        tt[2] = (t + t0) % 60;                  // Define SEC
        tt[0] = mmh / 60;                     // Define HOUR
        tt[1] = mmh % 60;                     // Define MIN
        for (int i = 0; i < 3; i++) {
            String s = Integer.toString(tt[i]);
            int len = s.length();
            if (len < 2) {
                s = "0" + tt[i];
            }
            if (i != 2) s += ":";
            tim.append(s);
        }
        return tim.toString();
    }


}
