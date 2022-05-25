"""GUI front-end to the Penetrance Estimator app."""

from pathlib import Path
from ast import literal_eval

from PyQt5 import uic, QtWidgets as QtW
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from cycler import cycler

from penEst import calc_stats, format_stats, isiter

#matplotlib.rcParams["text.usetex"] = True
        # turns out this was unnecessary, woo

# Canvas and MPL widget management inspired by/borrowed from
# https://stackoverflow.com/a/44029435

class PenEstCanvas(FigureCanvasQTAgg):
    """matplotlib canvas and figure for the Penetrance Estimator app."""
    
    def __init__(self):
        self.figure = Figure()
        self.ax_ftilde = self.figure.add_subplot(2, 1, 1)
        self.ax_ftildestar = self.figure.add_subplot(2, 1, 2)
        FigureCanvasQTAgg.__init__(self, self.figure)
        self.setSizePolicy(QtW.QSizePolicy.Expanding, QtW.QSizePolicy.Expanding)
        self.updateGeometry()
    
    def set_xaxis_var(self, varname):
        """Updates our plot to use a different X-axis variable."""
    
    def plot_data(self, stats):
        """Given the X-axis data and two lists of Y-axis datasets - one for
        ftilde, one for ftilde - clears our axes and plots the new data.
        
        """
        
        self.ax_ftilde.clear()
        self.ax_ftildestar.clear()
        
        # Which of our variables is an iterable (and thus our X-axis)?
        var_xaxis = "alpha" # default to alpha if none are iterable
        for varclass in ("alpha", "beta", "s"):
            var_xaxis = varclass if isiter(stats[varclass]) else var_xaxis
        
        self.restyleaxes(var_xaxis, stats[var_xaxis])
        
        # FIXME: change these in penEst.py eventually
        for yval in stats['stats']['betahat'].values():
            self.ax_ftilde.plot(stats[var_xaxis], yval)
        for yval in stats['stats']['betahatstar'].values():
            self.ax_ftildestar.plot(stats[var_xaxis], yval)
        
        self.draw()
    
    def restyleaxes(self, var_xaxis, xaxis_vals):
        """Sets up the axes styling our plots should have."""
        # Needed as separate routine because replotting resets them. Grr.
        
        # Axes names
        #self.ax_ftilde.set_xlabel(varname)
                # Not bothering to set this one as default layout covers it
                # anyways.
        self.ax_ftildestar.set_xlabel(var_xaxis)
        self.ax_ftilde.set_ylabel(r"$\tilde{f}$")
        self.ax_ftildestar.set_ylabel(r"$\tilde{f}^*$")
        
        # View limits
        # yaxis is always (0, 1). Usually the same is true of xaxis as well,
        # UNLESS we're using <s> as the xaxis. So we look at all current values
        # for x and take the largest one if that value > 1.
        bh_xmax = max(1, *xaxis_vals)
        self.ax_ftilde.set_xlim(xmin=0, xmax=bh_xmax)
        self.ax_ftildestar.set_xlim(xmin=0, xmax=bh_xmax)
        self.ax_ftilde.set_ylim(ymin=0, ymax=1)
        self.ax_ftildestar.set_ylim(ymin=0, ymax=1)
        
        # Property cyclers
        colors = cycler("color", ["k", "r", "g", "c", "gray", "m", "y", "b"])
        markers = cycler("marker", ["o", "x", "*", "s", "P", "D", "v", "^"])
        # We do the combination explicitly each time so that our plots don't
        # end up exhausting a single cycler. One cycler per plot!
        self.ax_ftilde.set_prop_cycle(colors + markers)
        self.ax_ftildestar.set_prop_cycle(colors + markers)

class PenEstPlot(QtW.QWidget):
    """Widget representing the space for a matplotlib canvas."""
    
    def __init__(self, parent):
        super().__init__(parent)
        # FIXME: at this point, remove the stylesheet
        self.canvas = PenEstCanvas()
        self.vbl = QtW.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)


UI_MainWindow = uic.loadUiType(str(
        Path(__file__).parent.joinpath("penestapp.ui")))[0]
class PenEstAppMain(QtW.QMainWindow, UI_MainWindow):
    """Main (only) window for a GUI program for calculating and plotting
    our penetrance estimations.
    
    """
    
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle("Penetrance Estimator")
        
        # Stats for the current plot
        self.stats = None
        
        # showing an initial plot (using values that are embedded in
        # penestapp.ui) so there's something shinyesque to look at
        self.refreshPlot()
    
    @property
    def newstats(self):
        """Returns the current vars and calculated stats for those vars."""
        
        modelid = self.VarModel.currentIndex() + 1  # because index starts at 0
        alpha = literal_eval(self.VarAlpha.text())
        beta = literal_eval(self.VarBeta.text())
        s = self.VarS.value()
        k = literal_eval(self.VarK.text())
        
        return calc_stats(modelid, k, alpha, beta, s)
    
    def refreshPlot(self):
        """Plots stats for the current entered variable values."""
        
        self.stats = self.newstats
        self.PlotWidget.canvas.plot_data(self.stats)
    
    def savePlotImage(self):
        """Saves the current plot to an image file."""
        
        filterstring = ' '.join([f"*.{ext}" for ext in
                self.PlotWidget.canvas.get_supported_filetypes().keys()])
        fileselect = QtW.QFileDialog.getSaveFileName(self,
                "Save Image", filter=f"Images ({filterstring})")
        
        if fileselect[0]:  #verify cancel wasn't pressed
            self.PlotWidget.canvas.figure.savefig(fileselect[0])
    
    def saveStatsFile(self):
        """Saves the current plotted data to a text file."""
        
        fileselect = QtW.QFileDialog.getSaveFileName(self,
                "Save Stats", filter="Text files (*.txt)")
        
        if fileselect[0]:
            with open(fileselect[0], 'w') as outfile:
                outfile.write(format_stats(self.stats))
