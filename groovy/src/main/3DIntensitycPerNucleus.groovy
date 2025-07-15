import ij.IJ
import ij.ImagePlus
import ij.gui.Roi
import ij.gui.ShapeRoi
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import ij.plugin.RGBStackMerge
import org.apache.commons.math3.stat.inference.TTest
import mcib3d.geom.Object3D
import mcib3d.geom.Objects3DPopulation
import mcib3d.image3d.ImageInt

// INPUT UI
//
//#@File(label = "Input Segmentation Files Directory", style = "directory") inputFilesSeg
//#@File(label = " Input Raw Files Directory", style = "directory") inputFilesRaw
//#@File(label = "output directory", style = "directory") outputDir
//#@Integer(label = "Nuclei channel", value = 2) nucleiChannel
//#@Integer(label = "Telomere channel", value = 0) telomereChannel
//#@Integer(label = "Marker channel", value = 1) markerChannel


// IDE
//
//
def inputFilesTrf1 = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/telomere")
def inputFilesNuclei = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/nuclei")
def inputFilesRaw = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/preprocess")
def outputDir = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output")
def nucleiChannel = 0.intValue()
def trf1Channel = 1.intValue()

//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -Input Seg Files Dir: " + inputFilesTrf1)
IJ.log("    -Input Raw Files Dir: " + inputFilesRaw)
IJ.log("    -output Dir: " + outputDir)
IJ.log("    -Nuclei Channel: " + nucleiChannel)
IJ.log("    -Telomere Channel: " + trf1Channel)
IJ.log("                                                           ");
/** Get files (images) from input directory */
def listOfFiles = inputFilesRaw.listFiles();
def wtRt = new ResultsTable()
def koRt = new ResultsTable()
def pValueRt = new ResultsTable()
def counterWt = 0.intValue()
def counterKo = 0.intValue()
//WT
def wtColocsN_All = new ArrayList<Double>()
/////Trf1
def wtTrf1N_All = new ArrayList<Double>()
def wtTrf1MeanInt_All = new ArrayList<Double>()
def wtTrf1SumInt_All = new ArrayList<Double>()
def wtTrf1StdInt_All = new ArrayList<Double>()

//KO
def koColocsN_All = new ArrayList<Double>()
//////Trf1
def koTrf1N_All = new ArrayList<Double>()
def koTrf1MeanInt_All = new ArrayList<Double>()
def koTrf1SumInt_All = new ArrayList<Double>()
def koTrf1StdInt_All = new ArrayList<Double>()
def wtTelomereIntensityMeanTotal = new ArrayList<Double>()
def wtTelomereIntensitySumTotal = new ArrayList<Double>()

for (def i = 0.intValue(); i < listOfFiles.length; i++) {
    if (listOfFiles[i].getName().contains("WT")) {

        def wtTelomereIntensityMeanPerNucleus = new ArrayList<Double>()
        def wtTelomereIntensitySumPerNucleus = new ArrayList<Double>()

        /** Get red labels */
        labelRed = new ImagePlus(inputFilesTrf1.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))

        /** Get nuclei image */
        def labelNuclei = new ImagePlus(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
        /** Get raw image */
        def imp = new ImagePlus(inputFilesRaw.getAbsolutePath() + File.separator + listOfFiles[i].getName())

        def channels = ChannelSplitter.split(imp)
        /** Get channels image */
        def chNuclei = channels[nucleiChannel.intValue()]
        def chRed = channels[trf1Channel.intValue()]

        // Get Nuclei objects population
        def imgNuclei = ImageInt.wrap(extractCurrentStack(labelNuclei));
        def populationNuclei = new Objects3DPopulation(imgNuclei);
        // Get Nuclei signal
        def signalNuclei = ImageInt.wrap(extractCurrentStack(chNuclei));

        // Get Red objects population
        def imgRed = ImageInt.wrap(extractCurrentStack(labelRed));
        def populationRed = new Objects3DPopulation(imgRed);
        // Get Red signal
        def signalRed = ImageInt.wrap(extractCurrentStack(chRed));


        for (int n = 0; n < populationNuclei.getNbObjects(); n++)
            for (int t = 0; t < populationRed.getNbObjects(); t++)
                if (populationNuclei.getObject(n).inside(populationRed.getObject(t).getCenterAsPoint())) {
                    wtTelomereIntensityMeanPerNucleus.add(populationRed.getObject(t).getPixMeanValue(signalRed))
                    wtTelomereIntensitySumPerNucleus.add(populationRed.getObject(t).getIntegratedDensity(signalRed))

                }
        wtTelomereIntensityMeanTotal.add(wtTelomereIntensityMeanPerNucleus.stream()
                .mapToDouble(d -> d)
                .average()
                .orElse(0.0))
        wtTelomereIntensitySumTotal.add(wtTelomereIntensitySumPerNucleus.stream()
                .mapToDouble(d -> d)
                .average()
                .orElse(0.0))
    }
}

def meanTelIntensityMean = wtTelomereIntensityMeanTotal.stream()
        .mapToDouble(d -> d)
        .average()
        .orElse(0.0)
def stdTelIntensityMean = std(wtTelomereIntensityMeanTotal, meanTelIntensityMean)
def meanTelIntensitySum = wtTelomereIntensitySumTotal.stream()
        .mapToDouble(d -> d)
        .average()
        .orElse(0.0)
def stdTelIntensitySum = std(wtTelomereIntensitySumTotal, meanTelIntensitySum)

def meanTelIntensityMeanp85 = meanTelIntensityMean + (1.04 * stdTelIntensityMean)
def meanTelIntensityMeanp15 = meanTelIntensityMean + (-1.04 * stdTelIntensityMean)

def meanTelIntensitySump85 = meanTelIntensitySum + (1.04 * stdTelIntensitySum)
def meanTelIntensitySump15 = meanTelIntensitySum + (-1.04 * stdTelIntensitySum)

for (def i = 0; i < listOfFiles.length; i++) {
    def counter = 0.intValue()
    def tablePerImage = new ResultsTable();
    /** Create image for each file in the input directory */
    def imps = new ImagePlus(inputFilesRaw.getAbsolutePath() + File.separator + listOfFiles[i].getName())
    def cal = imps.getCalibration()
    IJ.log(imps.getTitle())
    /** Split channels */
    def channels = ChannelSplitter.split(imps)

    /** Get trf1 channel */
    def labelTrf1 = new ImagePlus(inputFilesTrf1.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
    def chTrf1ToMeasure = channels[trf1Channel.intValue()]

    /** Get nuclei channel */
    def labelNuclei = new ImagePlus(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
    labelNuclei.setCalibration(cal)
    def chNucleiToMeasure = channels[nucleiChannel.intValue()]
    IJ.log(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))

    // Get Nuclei objects population
    def imgNuclei = ImageInt.wrap(extractCurrentStack(labelNuclei));
    def populationNuclei = new Objects3DPopulation(imgNuclei);
    // Get Nuclei signal
    def signalNuclei = ImageInt.wrap(extractCurrentStack(chNucleiToMeasure));

    // Get Trf1 objects population
    def imgTrf1 = ImageInt.wrap(extractCurrentStack(labelTrf1));
    def populationTrf1 = new Objects3DPopulation(imgTrf1);
    // Get Trf1 signal
    def signalTrf1 = ImageInt.wrap(extractCurrentStack(chTrf1ToMeasure));


    //IJ.saveAs(RGBStackMerge.mergeChannels(new ImagePlus[]{labelNuclei, chNucleiToMeasure, labelTrf1, chTrf1ToMeasure, labelH2ax, chH2axToMeasure}, false), "Tiff", outputDir.getAbsolutePath() + File.separator + "merge" + File.separator + listOfFiles[i].getName())

    if (listOfFiles[i].getName().contains("WT")) {
        //IJ.saveAs(RGBStackMerge.mergeChannels(new ImagePlus[]{labelNuclei, chNucleiToMeasure, labelTrf1, chTrf1ToMeasure, labelH2ax, chH2axToMeasure}, false), "Tiff", outputDir.getAbsolutePath() + File.separator + "merge" + File.separator + listOfFiles[i].getName())
        //WT clon Analysis
        def wtNucleiN = new ArrayList<Double>()
        //Area
        def wtNucleusArea = new ArrayList<Double>()
        //TRF1
        def wtTrf1N = new ArrayList<Double>()
        def wtTrf1MeanInt = new ArrayList<Double>()
        def wtTrf1SumInt = new ArrayList<Double>()
        def wtTrf1StdInt = new ArrayList<Double>()


        wtNucleiN.add(populationNuclei.getNbObjects().doubleValue())


        def nTelMorep85MeanTotal = new ArrayList<Double>()
        def nTelLowerp15MeanTotal = new ArrayList<Double>()
        def nTelMorep85SumTotal = new ArrayList<Double>()
        def nTelLowerp15SumTotal = new ArrayList<Double>()


        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
            def counterTrf1 = 0.intValue()

            def meanIntTrf1 = new ArrayList<Double>()
            def sumIntTrf1 = new ArrayList<Double>()
            def stdIntTrf1 = new ArrayList<Double>()
            def nTelMorep85Mean = new ArrayList<Object3D>()
            def nTelLowerp15Mean = new ArrayList<Object3D>()
            def nTelMorep85Sum = new ArrayList<Object3D>()
            def nTelLowerp15Sum = new ArrayList<Object3D>()

            wtNucleusArea.add(populationNuclei.getObject(j).volumeUnit)
            for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint())) {
                    counterTrf1++
                    meanIntTrf1.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    sumIntTrf1.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    stdIntTrf1.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                    //wtTrf1MeanInt_All.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    //wtTrf1SumInt_All.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    //wtTrf1StdInt_All.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                    if(populationTrf1.getObject(k).getPixMeanValue(signalTrf1) >= meanTelIntensityMeanp85)
                        nTelMorep85Mean.add(populationTrf1.getObject(k))
                    if(populationTrf1.getObject(k).getPixMeanValue(signalTrf1) <= meanTelIntensityMeanp15)
                        nTelLowerp15Mean.add(populationTrf1.getObject(k))
                    if(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1) >= meanTelIntensitySump85)
                        nTelMorep85Sum.add(populationTrf1.getObject(k))
                    if(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1) <= meanTelIntensitySump15)
                        nTelLowerp15Sum.add(populationTrf1.getObject(k))

                }
            }


            wtTrf1N.add(counterTrf1.doubleValue())
            wtTrf1N_All.add(counterTrf1.doubleValue())
            wtTrf1MeanInt.add(meanIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtTrf1MeanInt_All.add(meanIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtTrf1SumInt.add(sumIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .sum())
            wtTrf1SumInt_All.add(sumIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtTrf1StdInt.add(stdIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtTrf1StdInt_All.add(stdIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))

            nTelMorep85MeanTotal.add(nTelMorep85Mean.size())
            nTelLowerp15MeanTotal.add(nTelLowerp15Mean.size())
            nTelMorep85SumTotal.add(nTelMorep85Sum.size())
            nTelLowerp15SumTotal.add(nTelLowerp15Sum.size())


        }
        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {

            def zRange = (populationNuclei.getObject(j).zmax - populationNuclei.getObject(j).zmin)
            if (zRange > 2) {
                counterWt++
                wtRt.incrementCounter()
                wtRt.setValue("Image Title", counterWt, listOfFiles[i].getName())
                wtRt.setValue("N of Nuclei per Image", counterWt, populationNuclei.getNbObjects())
                wtRt.setValue("Nucleus Label", counterWt, populationNuclei.getObject(j).getValue())
                wtRt.setValue("Nucleus Volume (microns3)", counterWt, populationNuclei.getObject(j).getVolumeUnit())
                wtRt.setValue("Mean of Nucleus Volume (microns3) per Image", counterWt, wtNucleusArea.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))
                wtRt.setValue("N of Telomeres per Nucleus", counterWt, wtTrf1N.get(j).doubleValue())
                wtRt.setValue("N of Telomeres per Nucleus >p85mean", counterWt, nTelMorep85MeanTotal.get(j).doubleValue())
                wtRt.setValue("N of Telomeres per Nucleus <p15mean", counterWt, nTelLowerp15MeanTotal.get(j).doubleValue())
                wtRt.setValue("N of Telomeres per Nucleus >p85sum", counterWt, nTelMorep85SumTotal.get(j).doubleValue())
                wtRt.setValue("N of Telomeres per Nucleus <p15sum", counterWt, nTelLowerp15SumTotal.get(j).doubleValue())
                wtRt.setValue("p15 Mean Intensity", counterWt, meanTelIntensityMeanp15)
                wtRt.setValue("p85 Mean Intensity", counterWt,meanTelIntensityMeanp85)
                wtRt.setValue("p15 Sum Intensity", counterWt, meanTelIntensitySump15)
                wtRt.setValue("p85 Sum Intensity", counterWt, meanTelIntensitySump85)


                wtRt.setValue("Mean Intensity of Telomeres per Nucleus", counterWt, wtTrf1MeanInt.get(j))

                wtRt.setValue("Sum Intensity of Telomeres per Nucleus", counterWt, wtTrf1SumInt.get(j))


                wtRt.setValue("Mean of Sum Intensity of Telomeres per Nucleus", counterWt, wtTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))
                wtRt.setValue("Std of Sum Intensity of Telomeres per Nucleus", counterWt, std(wtTrf1SumInt, wtTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0)))
                wtRt.setValue("SEM of Sum Intensity of Telomeres per Nucleus", counterWt, std(wtTrf1SumInt, wtTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0)) / Math.sqrt(wtTrf1SumInt.size()))
            }

        }

    }


    if (listOfFiles[i].getName().contains("KO")) {
        //IJ.saveAs(RGBStackMerge.mergeChannels(new ImagePlus[]{labelNuclei, chNucleiToMeasure, labelTrf1, chTrf1ToMeasure, labelH2ax, chH2axToMeasure}, false), "Tiff", outputDir.getAbsolutePath() + File.separator + "merge" + File.separator + listOfFiles[i].getName())
        //KO clon Analysis
        def koNucleiN = new ArrayList<Double>()
        //Area
        def koNucleusArea = new ArrayList<Double>()
        //TRF1
        def koTrf1N = new ArrayList<Double>()
        def koTrf1MeanInt = new ArrayList<Double>()
        def koTrf1SumInt = new ArrayList<Double>()
        def koTrf1StdInt = new ArrayList<Double>()


        koNucleiN.add(populationNuclei.getNbObjects().doubleValue())
        def nTelMorep85MeanTotal = new ArrayList<Double>()
        def nTelLowerp15MeanTotal = new ArrayList<Double>()
        def nTelMorep85SumTotal = new ArrayList<Double>()
        def nTelLowerp15SumTotal = new ArrayList<Double>()

        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
            def counterTrf1 = 0.intValue()
            def meanIntTrf1 = new ArrayList<Double>()
            def sumIntTrf1 = new ArrayList<Double>()
            def stdIntTrf1 = new ArrayList<Double>()

            def nTelMorep85Mean = new ArrayList<Object3D>()
            def nTelLowerp15Mean = new ArrayList<Object3D>()
            def nTelMorep85Sum = new ArrayList<Object3D>()
            def nTelLowerp15Sum = new ArrayList<Object3D>()

            koNucleusArea.add(populationNuclei.getObject(j).volumeUnit)
            for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint())) {
                    counterTrf1++
                    meanIntTrf1.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    sumIntTrf1.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    stdIntTrf1.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                    //koTrf1MeanInt_All.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    //koTrf1SumInt_All.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    //koTrf1StdInt_All.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                    if(populationTrf1.getObject(k).getPixMeanValue(signalTrf1) >= meanTelIntensityMeanp85)
                        nTelMorep85Mean.add(populationTrf1.getObject(k))
                    if(populationTrf1.getObject(k).getPixMeanValue(signalTrf1) <= meanTelIntensityMeanp15)
                        nTelLowerp15Mean.add(populationTrf1.getObject(k))
                    if(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1) >= meanTelIntensitySump85)
                        nTelMorep85Sum.add(populationTrf1.getObject(k))
                    if(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1) <= meanTelIntensitySump15)
                        nTelLowerp15Sum.add(populationTrf1.getObject(k))
                }
            }


            koTrf1N.add(counterTrf1.doubleValue())
            koTrf1N_All.add(counterTrf1.doubleValue())
            koTrf1MeanInt.add(meanIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koTrf1MeanInt_All.add(meanIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koTrf1SumInt.add(sumIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .sum())
            koTrf1SumInt_All.add(sumIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .sum())
            koTrf1StdInt.add(stdIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koTrf1StdInt_All.add(stdIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            nTelMorep85MeanTotal.add(nTelMorep85Mean.size())
            nTelLowerp15MeanTotal.add(nTelLowerp15Mean.size())
            nTelMorep85SumTotal.add(nTelMorep85Sum.size())
            nTelLowerp15SumTotal.add(nTelLowerp15Sum.size())

        }
        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {


            def zRange = (populationNuclei.getObject(j).zmax - populationNuclei.getObject(j).zmin)
            if (zRange > 2) {
                counterKo++
                koRt.incrementCounter()
                koRt.setValue("Image Title", counterKo, listOfFiles[i].getName())
                koRt.setValue("Nucleus Label", counterKo, populationNuclei.getObject(j).getValue())
                koRt.setValue("N of Nuclei per Image", counterKo, populationNuclei.getNbObjects())
                koRt.setValue("Nucleus Volume (microns3)", counterKo, populationNuclei.getObject(j).getVolumeUnit())
                koRt.setValue("Mean of Nucleus Volume (microns3) per Image", counterKo, koNucleusArea.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))
                IJ.log(counterKo + "-----")
                koRt.setValue("N of Telomeres per Nucleus", counterKo, koTrf1N.get(j).doubleValue())
                koRt.setValue("N of Telomeres per Nucleus >p85mean", counterKo, nTelMorep85MeanTotal.get(j).doubleValue())
                koRt.setValue("N of Telomeres per Nucleus <p15mean", counterKo, nTelLowerp15MeanTotal.get(j).doubleValue())
                koRt.setValue("N of Telomeres per Nucleus >p85sum", counterKo, nTelMorep85SumTotal.get(j).doubleValue())
                koRt.setValue("N of Telomeres per Nucleus <p15sum", counterKo, nTelLowerp15SumTotal.get(j).doubleValue())
                koRt.setValue("p15 Mean Intensity", counterKo, meanTelIntensityMeanp15)
                koRt.setValue("p85 Mean Intensity", counterKo,meanTelIntensityMeanp85)
                koRt.setValue("p15 Sum Intensity", counterKo, meanTelIntensitySump15)
                koRt.setValue("p85 Sum Intensity", counterKo, meanTelIntensitySump85)
                koRt.setValue("Mean Intensity of Telomeres per Nucleus", counterKo, koTrf1MeanInt.get(j))
                koRt.setValue("Sum Intensity of Telomeres per Nucleus", counterKo, koTrf1SumInt.get(j))
                koRt.setValue("Mean of Sum Intensity of Telomeres per Nucleus", counterKo, koTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))
                koRt.setValue("Std of Sum Intensity of Telomeres per Nucleus", counterKo, std(koTrf1SumInt, koTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0)))
                koRt.setValue("SEM of Sum Intensity of Telomeres per Nucleus", counterKo, std(koTrf1SumInt, koTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0)) / Math.sqrt(koTrf1SumInt.size()))

            }

        }
    }


}

wtRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_WT_per_Nucleus_percentile.csv")
koRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_KO_per_Nucleus_percentile.csv")

ImagePlus extractCurrentStack(ImagePlus plus) {
    // check dimensions
    int[] dims = plus.getDimensions();//XYCZT
    int channel = plus.getChannel();
    int frame = plus.getFrame();
    ImagePlus stack;
    // crop actual frame
    if ((dims[2] > 1) || (dims[4] > 1)) {
        IJ.log("hyperstack found, extracting current channel " + channel + " and frame " + frame);
        def duplicator = new Duplicator();
        stack = duplicator.run(plus, channel, channel, 1, dims[3], frame, frame);
    } else stack = plus.duplicate();

    return stack;
}

static double std(ArrayList<Double> table, double mean) {
    // Step 1:
    double meanDef = mean
    double temp = 0;

    for (int i = 0; i < table.size(); i++) {
        int val = table.get(i);

        // Step 2:
        double squrDiffToMean = Math.pow(val - meanDef, 2);

        // Step 3:
        temp += squrDiffToMean;
    }

    // Step 4:
    double meanOfDiffs = (double) temp / (double) (table.size());

    // Step 5:
    return Math.sqrt(meanOfDiffs);
}