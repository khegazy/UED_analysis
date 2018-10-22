import sys
sys.argv.append("-b")
import ROOT


class PLOTclass :

  # Varaible Members
  genStyle = ROOT.TStyle("genStyle", "genHistStyle")

  #Function Members
  @staticmethod
  def centerAxisTitles(hist) :
    hist.GetXaxis().CenterTitle(1)
    hist.GetYaxis().CenterTitle(1)

    return




  # Constructor
  def __init__(self) :

    # use plain black on white colors
    icol=0; # WHITE
    self.genStyle.SetFrameBorderMode(icol);
    self.genStyle.SetFrameFillColor(icol);
    self.genStyle.SetCanvasBorderMode(icol);
    self.genStyle.SetCanvasColor(icol);
    self.genStyle.SetPadBorderMode(icol);
    self.genStyle.SetPadColor(icol);
    self.genStyle.SetStatColor(icol);

    # set the paper & margin sizes
    self.genStyle.SetPaperSize(20,26);

    # set margin sizes
    self.genStyle.SetPadTopMargin(0.04);
    #self.genStyle.SetPadRightMargin(0.05);
    self.genStyle.SetPadRightMargin(0.15);
    self.genStyle.SetPadBottomMargin(0.14);
    self.genStyle.SetPadLeftMargin(0.14);

    # set title offsets (for axis label)
    self.genStyle.SetTitleX(0.5)
    #self.genStyle.SetTitleY(0.5)
    self.genStyle.SetTitleW(0.6)
    self.genStyle.SetTitleAlign(23)
    #self.genStyle.SetTitleXOffset(1.4);
    #self.genStyle.SetTitleYOffset(1.4);

    # use large fonts
    # Int_t font=72; // Helvetica italics

    font=42; # Helvetica
    tsize=0.04;
    self.genStyle.SetTextFont(font);

    self.genStyle.SetTextSize(tsize);
    self.genStyle.SetLabelFont(font,"x");
    self.genStyle.SetTitleFont(font,"x");
    self.genStyle.SetLabelFont(font,"y");
    self.genStyle.SetTitleFont(font,"y");
    self.genStyle.SetLabelFont(font,"z");
    self.genStyle.SetTitleFont(font,"z");

    self.genStyle.SetLabelSize(tsize,"x");
    self.genStyle.SetTitleSize(tsize,"x");
    self.genStyle.SetLabelSize(tsize,"y");
    self.genStyle.SetTitleSize(tsize,"y");
    self.genStyle.SetLabelSize(tsize,"z");
    self.genStyle.SetTitleSize(tsize,"z");

    # use bold lines and markers
    self.genStyle.SetMarkerStyle(20);
    self.genStyle.SetMarkerSize(1.2);
    self.genStyle.SetHistLineWidth(2);
    self.genStyle.SetLineStyleString(2,"[12 12]");  #postscript dashes

    # get rid of X error bars 
    #self.genStyle.SetErrorX(0.001);
    # get rid of error bar caps
    self.genStyle.SetEndErrorSize(0.);

    # do not display any of the standard histogram decorations
    self.genStyle.SetOptTitle(0);
    #self.genStyle.SetOptStat(1111);
    self.genStyle.SetOptStat(0);
    #self.genStyle.SetOptFit(1111);
    self.genStyle.SetOptFit(0);

    # put tick marks on top and RHS of plots
    self.genStyle.SetPadTickX(1);
    self.genStyle.SetPadTickY(1);

    # set color palette 
    self.genStyle.SetPalette(55);
    #self.genStyle.SetPalette(kDarkBodyRadiator);
    #self.genStyle.SetPalette(kVisibleSpectrum);

    return

