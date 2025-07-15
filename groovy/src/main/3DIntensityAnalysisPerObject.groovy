import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import ij.plugin.RGBStackMerge
import inra.ijpb.algo.AlgoListener
import inra.ijpb.label.LabelImages
import inra.ijpb.measure.ResultsBuilder
import inra.ijpb.morphology.Morphology
import inra.ijpb.morphology.Strel3D
import inra.ijpb.morphology.strel.CubeStrel
import inra.ijpb.plugins.AnalyzeRegions
import inra.ijpb.plugins.AnalyzeRegions3D
import inra.ijpb.measure.IntensityMeasures
import mcib3d.geom.Object3D
import mcib3d.geom.Objects3DPopulation
import mcib3d.image3d.ImageInt


import java.io.File


// INPUT UI
//
//#@File(label = "Input File Directory Telomere", style = "directory") inputFilesNuclei
//#
//@File(label = "Input Raw Files Directory", style = "directory") inputFiles
//#
//@File(label = "Output directory", style = "directory") outputDir


// IDE
//
//
//def headless = true;
//new ij.ImageJ().setVisible(false);


def inputFilesRed = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_14_jfiel/output/telomere")
def inputFilesNuclei = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_14_jfiel/output/nuclei")
def inputFilesRaw = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_14_jfiel/output/images")
def outputDir = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_14_jfiel/output/")
def nucleiChannel = 0.intValue()
def trf1Channel = 1.intValue()

IJ.log("-Parameters selected: ")
IJ.log("    -inputFilesNuclei: " + inputFilesNuclei.toString())
IJ.log("    -inputFilesRedGreen: " + inputFilesRed.toString())
IJ.log("    -inputRawFiles: " + inputFilesRaw.toString())
IJ.log("    -outputDir: " + outputDir.toString())
IJ.log("                                                           ");


def listOfFiles = inputFilesRaw.listFiles();
int row = 0;

def wtRt = new ResultsTable()
def koRt = new ResultsTable()

def counterWt = 0.intValue()
def counterKo = 0.intValue()


for (def i = 0.intValue(); i < listOfFiles.length; i++) {
    if (listOfFiles[i].getName().contains("WT")) {
        /** Get red labels */
        labelRed = new ImagePlus(inputFilesRed.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))

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


        for (int n = 0; n < populationNuclei.getNbObjects(); n++) {


            for (int t = 0; t < populationRed.getNbObjects(); t++) {
                if (populationNuclei.getObject(n).inside(populationRed.getObject(t).getCenterAsPoint())) {
                    counterWt++
                    wtRt.incrementCounter();
                    int[] numbersRed = populationNuclei.getObject(n).getNumbering(imgRed)
                    wtRt.setValue("Image", counterWt, listOfFiles[i].getName())
                    wtRt.setValue("Telomere ID", counterWt, populationRed.getObject(t).getValue())
                    wtRt.setValue("Telomere Mean Intensity per Telomere", counterWt, populationRed.getObject(t).getPixMeanValue(signalRed))
                    wtRt.setValue("Telomere Sum Intensity per Telomere", counterWt, populationRed.getObject(t).getIntegratedDensity(signalRed))
                    wtRt.setValue("Telomere Std Intensity per Telomere", counterWt, populationRed.getObject(t).getPixStdDevValue(signalRed))
                    wtRt.setValue("Telomere Volume per Telomere", counterWt, populationRed.getObject(t).getVolumeUnit())
                    wtRt.setValue("Nucleus ID Label ", counterWt, populationNuclei.getObject(n).getValue())
                    //wtRt.setValue("N of Telomere per Nucleus", counterWt, numbersRed[0])
                    //wtRt.setValue("Telomere Volume occupied per Nucleus", counterWt, numbersRed[1])

                }
            }


        }
    }

    if (listOfFiles[i].getName().contains("KO")) {
        /** Get red labels */
        labelRed = new ImagePlus(inputFilesRed.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))

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


        for (int n = 0; n < populationNuclei.getNbObjects(); n++) {


            for (int t = 0; t < populationRed.getNbObjects(); t++) {
                if (populationNuclei.getObject(n).inside(populationRed.getObject(t).getCenterAsPoint())) {
                    counterKo++
                    koRt.incrementCounter();
                    int[] numbersRed = populationNuclei.getObject(n).getNumbering(imgRed)
                    koRt.setValue("Image", row, listOfFiles[i].getName())
                    koRt.setValue("Telomere ID", row, populationRed.getObject(t).getValue())
                    koRt.setValue("Telomere Mean Intensity per Telomere", row, populationRed.getObject(t).getPixMeanValue(signalRed))
                    koRt.setValue("Telomere Sum Intensity per Telomere", row, populationRed.getObject(t).getIntegratedDensity(signalRed))
                    koRt.setValue("Telomere Std Intensity per Telomere", row, populationRed.getObject(t).getPixStdDevValue(signalRed))
                    koRt.setValue("Telomere Volume per Telomere", row, populationRed.getObject(t).getVolumeUnit())
                    koRt.setValue("Nucleus ID Label ", row, populationNuclei.getObject(n).getValue())
                    //koRt.setValue("N of Telomere per Nucleus", row, numbersRed[0])
                    //koRt.setValue("Telomere Volume occupied per Nucleus", row, numbersRed[1])

                }
            }


        }
    }

}
//Save results table per cell or volume as csv file
wtRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_WT_per_Telomere.csv")
koRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_KO_per_Telomere.csv")

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

static double std(ArrayList<Double> table, double meanArea) {
    // Step 1:
    double mean = meanArea;
    double temp = 0;

    for (int i = 0; i < table.size(); i++) {
        double val = table.get(i);

        // Step 2:
        double squrDiffToMean = Math.pow(val - mean, 2);

        // Step 3:
        temp += squrDiffToMean;
    }

    // Step 4:
    double meanOfDiffs = (double) temp / (double) (table.size());

    // Step 5:
    return Math.sqrt(meanOfDiffs);
}

IJ.log("Done!!!")



