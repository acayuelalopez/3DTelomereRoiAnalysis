import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.WindowManager
import ij.measure.Calibration
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import ij.plugin.frame.RoiManager
import inra.ijpb.label.LabelImages
import inra.ijpb.measure.region3d.RegionAnalyzer3D
import inra.ijpb.morphology.Strel
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import mcib3d.geom.Objects3DPopulation
import mcib3d.image3d.ImageInt
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest
import java.util.stream.Stream;

// INPUT UI
//
//#@File(label = "Input LIF File Directory", style = "directory") inputFiles
//#@File(label = "Output directory", style = "directory") outputDir
//#@File(label = "Green model", style = "file") greenModel
//#@File(label = "Red model", style = "file") redModel
//#@Integer(label = "Reference Channel", value = 1) refIndex
//#@Integer(label = "Target Channel", value = 2) targetIndex
def inputFiles = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/images")
def outputDir = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/preprocess")
def rm = new RoiManager()
rm = RoiManager.getInstance()
// IDE
//
//
//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -inputFileDir: " + inputFiles)
IJ.log("    -outputDir: " + outputDir)

IJ.log("                                                           ");
/** Get files (images) from input directory */
def listOfFiles = inputFiles.listFiles(); ;


def rois = rm.getRoisAsArray()
for (def k = 0; k < rois.length; k++) {
    IJ.log(rm.getName(rm.getRoiIndex(rois[k])).replaceAll("/", ""))
    def roiName = rm.getName(rm.getRoiIndex(rois[k])).replaceAll("/", "")
    def imp = new ImagePlus(inputFiles.getAbsolutePath() + File.separator + rm.getName(rm.getRoiIndex(rois[k])).replaceAll("/", "") + ".tif")
    // Create a new stack to hold the processed slices
    def newStack = new ImageStack(imp.getWidth(), imp.getHeight())
    // Iterate through slices of the z-stack
    def stack = imp.getStack()
    for (int i = 1; i <= stack.getSize(); i++) {
        // Get the processed slice as an ImagePlus
        def slice = new ImagePlus("slice" + i, stack.getProcessor(i))
        slice.setRoi(rois[k])
        IJ.run(slice, "Make Inverse", "")
        IJ.setBackgroundColor(0, 0, 0)
        IJ.run(slice, "Clear", "slice")
        newStack.addSlice(slice.getProcessor())
    }
    IJ.saveAs(imp, "Tiff", outputDir.getAbsolutePath() + File.separator + rm.getName(rm.getRoiIndex(rois[k])).replaceAll("/", ""))


}
