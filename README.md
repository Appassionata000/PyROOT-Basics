# PyROOT Basics

This tutorial is based on the six-week undergraduate Third Year lab work performed by at the University of Manchester. The tutorial focuses on data fitting and plot configuration and is a summary of codes useful during the lab. The root files that pass the selection cuts are obtained using codes written by senior students. The repository owner does not have permission to upload these files. The original data are from [ATLAS Open Data 13 TeV](http://opendata.atlas.cern/release/2020/documentation/index.html). 

**Installing ROOT on MacOS**

```shell
# at terminal
brew install ROOT
vim ~/.bash_profile
# or vim ~/.zshrc depending on your shell
source /usr/local/bin/thisroot.sh  # add in .bash_profile or .zshrc
```

```python
# Import ROOT in Python
import ROOT as r
```

## Reading files and plotting histograms

```python
# create a canvas.
canvas1 = r.TCanvas("histograms_1")		
f_photonPt_Cuts1 = r.TFile.Open("out/ATLAS_yy_Cuts1.root", "READ") # hfile1 is a TFile object
f_photonPt_Cuts2 = r.TFile.Open("out/ATLAS_yy_Cuts2.root", "READ")
f_photonPt_Cuts3 = r.TFile.Open("out/ATLAS_yy_Cuts3.root", "READ")
# Can use print(f_photonPt_cuts1.ls()) to list the keys in the file
# Get histogram data using their keys
photonPt_Cuts1 = f_photonPt_Cuts1.Get("photon_pt[0]")	# photonPt_1 is a TH1F object
photonPt_Cuts2 = f_photonPt_Cuts2.Get("photon_pt[0]")
photonPt_Cuts3 = f_photonPt_Cuts3.Get("photon_pt[0]")
# Avoid losing data after closing the file
photonPt_1.SetDirectory(0) 
photonPt_2.SetDirectory(0)
photonPt_3.SetDirectory(0)
f_photonPt_Cuts1.Close()	
f_photonPt_Cuts2.Close()	
f_photonPt_Cuts3.Close()	
... # Analysis/Fitting/...
photonPt_Cuts1.Draw()
photonPt_Cuts2.Draw("SAME")  # Draw in the same canvas
photonPt_Cuts3.Draw("SAME")
```

## Fitting

The `TH1::Fit()` method:

```cpp
// Fit a TF1 function with the TH1::Fit() method
TFitResultPtr Fit(TF1 *f1, const char *option="",    // commonly used options are listed below
						   const char *goption="",					   // the same options as TF1::Draw()
						double xxmin=0.0, double xxmax=0.0);   // Specify the range over which to apply the fit
// Fit a pre-defined function with the TH1::Fit() method
TFitResultPtr Fit(const char *fname, const char *option="", const char *goption="", 
						double xxmin=0.0, double xxmax=0.0);
// Get list of pre-defined functions (fname)
TF1::InitStandardFunctions();
gROOT->GetListOfFunctions()->ls();
// Predefined functions: "gaus", "gausn", "expo", "pol N", "chebyshev", "landau".
```

Using `Fit()` in PyROOT with a pre-defined function can simply be done by a single line:

```python
invMass.Fit("gaus", "Q") # Gaussian fit with no information displayed after fitting
```

Commonly used fit options

| `option` | Effect                                                       |
| -------- | ------------------------------------------------------------ |
| `"L"`    | Use log likelihood method. To be used when the histogram represents counts. |
| `"Q"`    | Quiet mode (minimum printing)                                |
| `"E"`    | Performing better errors                                     |
| `"S"`    | The result of the fit is returned in the `TFitResultPtr`, otherwise a `nullptr` |
| `"M"`    | Improve fit results by using the improve algorithm of  `TMinuit` |
| `"R"`    | Use the range specified in the function range                |
| `"N"`    | Do not store the graphics function, do not draw              |
| `"0"`    | Do not plot the result of the fit                            |
| `"+"`    | Add this new fitted function to the list of fitted functions (By default the previous  function is deleted) |
| `"B"`    | Use this option when you want to fix one or more parameters and the fitting function is a predefined one |
| `"C"`    | In case of linear fitting, don't calculate the chisquare (saves time) |
| `"F"`    | If fitting a linear function, switch to use the default minimizer. By default, `polN` are fitted by the linear fitter |

### User-defined fitting functions

One-dimensional functions in ROOT are typically defined as `TF1` objects

```cpp
// A common constructor and examples of initializing TF1 objects (C++)
TF1::TF1(const char *name, const char *formula, Double_t xmin = 0, Double_t xmax = 1, Option_t *option)	
auto func1 = new TF1("func1", "sin(x)/x", 0, 10);
auto func2 = new TF1("func2", "TMath::BesselJ0(x)", 0, 10);
auto quartic = new TF1("quartic", "[0]*x**4 + [1]*x**3 + [2]*x**2 + [3]*x + [4]", 0, 10);
```

```python
# Corresponding TF1 objects in Python
func1 = r.TF1("func1", "sin(x)/x", 0, 10)
# The 2nd parameter of the TF1 constructor is a text string parsed as a TFormla. It should obey C++ syntax instead of using TMath.BesselJ0(x)
func2 = r.TF1("func2", "TMath::BesselJ0(x)", 0, 10)
quartic = r.TF1("quartic", "[0]*x**4 + [1]*x**3 + [2]*x**2 + [3]*x + [4]")
```

- Pre-defined functions can be found in the `TFormula` class.

#### Define more complicated functions 

```python
# The function always takes two parameters no matter how many fitting parameters there are
# and the variable x should be written as x[0]
def quartic_tot(x, p):
    return p[0] * x[0]**4 + p[1] * x[0]**3 + p[2] * x[1]**2 + p[3] * x[0] + p[4]
```

```python
# A function with reject region can used to fit the background sidebands
def quartic_bkg(x, p):    
    if (x[0] >= 120e3 and x[0] <= 130e3):
       r.TF1.RejectPoint(True)
    return p[0] + p[1] * x[0] + p[2] * x[0]**2 + p[3] * x[0]**3 + p[4] * x[0]**4
```

```python
import math
def double_sided_crystal_ball(x, p):
	# Double sided crystal ball function can be used to fit the signal accurately
    # p[0]: constant, p[1]: mean, p[2]: std,   p[3]: a1(low),  p[4]: n1(low),  p[5]: a2(high)  p[6]: n2(high)
    s = (x[0] - p[1]) / p[2]
    a1, n1, a2, n2 = p[3], p[4], p[5], p[6]
    result = 1     
    if (n1 / abs(a1) - abs(a1) - s) > 0 and (n2 / abs(a2) - abs(a2) + s) > 0:
    # This if statement is essential, since ROOT fitting process might come across 
    # undefined math value, such as the square root of a negative number.
        A1 = math.pow((n1 / abs(a1)), n1) * math.exp(-a1**2 / 2) 
        A2 = math.pow((n2 / abs(a2)), n2) * math.exp(-a2**2 / 2)
        B1 = math.pow((n1 / abs(a1) - abs(a1) - s), -n1)
        B2 = math.pow((n2 / abs(a2) - abs(a2) + s), -n2)
        if s < -a1:
            result = p[0] * A1 * B1
        elif s > a2:
            result = p[0] * A2 * B2
        else:
            result = p[0] * math.exp(-s**2 / 2)
    return result
```


Standard process of fitting `TH1` object data with a user defined function in PyROOT

```python
# Example of fitting the MC signal using the double sided crystal ball function
# If the user wants to fit multiple functions to a set of data, it's a good idea to copy the original TH1 object first. To begin, if more fits are to be performed later, the first fitted result will be deleted by default. Second, when plotting, the fitting result will always be on top of the TH1 data.
invMass_0 = invMass.Clone("Original data")
invMass_dscb = invMass.Clone("DSCB fit")
# Define the fitting function TF1 object and specify the fitting range as well as the parameters
dscb = r.TF1("dscb", double_sided_crystal_ball, 105e3, 140e3, 7)
# Essential! Setting initial parameters
dscb.SetParameters(1, 125e3, 2e3, 1, 10, 1, 20) 
# Setting boundaries for parameters also has an influence on ROOT fitting performance, especially for complicated functions like DSCB
dscb.SetParLimits(0, 0, 100)
dscb.SetParLimits(3, 0.1, 2)  # a1
dscb.SetParLimits(4, 1, 30)   # n1
dscb.SetParLimits(5, 0.1, 2)  # a2
dscb.SetParLimits(6, 1, 40)   # n2
# Normalize the data (optional)
norm = invMass.Integral()
invMass_0.Scale(1 / norm)
invMass_dscb.Scale(1 / norm)
# Perform the fit and get the fitted results
fit_result = invMass_dscb.Fit("dscb", "M, R, S, Q")
fitted_dscb = invMass_dscb.GetFunction("dscb")  # Fitted function
# Store fitted results
dscb_fit_pars = fitted_dscb.GetParameters()
dscb_fit_par_errors = fitted_dscb.GetParErrors()
chisquared = fit_result.MinFcnValue()
ndf = fit_result.Ndf()
# Get integral of the fitted function
int_func_sig = fitted_dscb.Integral(120e3, 130e3)
# Plotting 
invMass_0.Draw("E") # Draw with errorbars
fitted_dscb.Draw("SAME")
canvas.Draw()
canvas.SaveAs("fit.pdf")
```

## Configuring plots

- Plotstyle-related classes: [TAttMarker](https://root.cern.ch/doc/master/classTAttMarker.html), [TattLine](https://root.cern.ch/doc/master/classTAttLine.html), [TAttAxis](https://root.cern.ch/doc/master/classTAttAxis.html)

```python
# Global setting
r.gStyle.SetErrorX(0)
r.gStyle.SetOptStat(0) # don't display statistics box
# Canvas size
r.gStyle.SetCanvasDefH(height) # default 500
r.gStyle.SetCanvasDefW(width) # default 700
```
```python
# Pad
r.gPad.SetMargin(0.1, 0.2, 0.2, 0.1)  # left, right, bottom, top
r.gPad.SetFrameLineWidth(3) # Frame line width
r.gPad.SetLogy(1) # logrithm y-axis
```
```python
# TH1F/TF1/... objects setting
invMass.SetTitle("Invariant Mass for the diphoton system")
```
```python
# Axis label
invMass.SetLabelSize(0.07, "X") # Relative size
invMass.SetLabelSize(0.07, "Y")
invMass.GetXaxis().SetNdivisions(10) # Set the division of X axis
invMass.GetYaxis().SetNdivisions(10)
invMass.GetXaxis().SetTitle("Invariant Mass") # Set label titles
invMass.GetYaxis().SetTitle("Entries")
invMass.GetXaxis().SetTitleSize(0.06)
invMass.GetYaxis().SetTitleSize(0.06)
invMass.SetTitleOffset(1, "X") # Modify positions
invMass.SetTitleOffset(1, "Y")
```
```python
# Marker
invMass.SetMarkerStyle(8)
```
```python
# Line
fitted_bkg.SetLineStyle(2)
fitted_bkg.SetLineColorAlpha(4, 0.9)
fitted_bkg.SetLineWidth(3)
```
```python
# Axis range
invMass.GetYaxis().SetRangeUser(0, 8000)
```
```python
# Set statistics box
r.gStyle.SetOptStat(1111)
```

```python
# Set statistics box position
var.Draw() # Set statistics box position after Draw()
r.gPad.Update()
pave = var.FindObject("stats")
pave.SetTextSize(0.055)
pave.SetX1NDC(0.8)
pave.SetY1NDC(0.2)
pave.SetX2NDC(0.95)
pave.SetY2NDC(0.9)
```

```python
# Create a TLegend object
legend = r.TLegend(0.55, 0.55,	# relative coordinate (x1, y1) (left down corner)
                   0.85, 0.85)  # relative coordinate (x2, y2) (right up corner)
legend.AddEntry(ATLAS_DATA, " Data")
legend.AddEntry(MC_DATA, " Signal")
legend.AddEntry(BackgroundFit, " Background fit")
legend.AddEntry(InclusiveFit, " Signal + Background fit")
legend.SetBorderSize(0)
r.gStyle.SetLegendTextSize(0.04)
legend.Draw("SAME")
```

```python
# Create a TLatex object
infolatex1 = r.TLatex()
infolatex1.SetTextSize(0.036)	# Relative font size 
infolatex1.DrawLatexNDC(0.15, 0.36,  "\sqrt{s} = 13 TeV, \int Ldt=10.064 fb^{-1}") # Relative position
infolatex1.Draw("SAME")
```

## Miscellaneous

```python
# Get list of errors
original_errors = []
for i in range(invMass.GetNbinsX()):
    original_errors = original_errors + [invMass.GetBinError(i + 1)]
```

```python
# Subtract the data from the background fitting
for i in range(invMass.GetNbinsX()):
    invMass_sub.AddBinContent(i+1, -1 * (fitted_bkg(invMass_sub.GetBinCenter(i+1))))
    #invMass_sub.SetBinError(i+1, original_errors[i]) # reset errors
```

```python
# TF1 integral
bin_width = 2e3
N_fit_sig = fitted_gaus.Integral(120e3, 130e3) / bin_width
N_fit_sig_err = fitted_gaus.IntegralError(120e3, 130e3) / bin_width
```

```python
# TH1F integral
import ctypes
N_selec = invMass.Integral(11, 15) # bin number
N_selec_err = c_double()
invMass.IntegralAndError(11, 15, N_selec_err)
```

































































