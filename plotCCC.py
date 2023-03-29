#!/usr/bin/env python

import CombineHarvester.CombineTools.plotting as plot
import ROOT
import argparse
import json

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()
ROOT.gStyle.SetTickLength(0., "Y")

parser = argparse.ArgumentParser()
parser.add_argument('--input','-i',help = 'input json file')
parser.add_argument('--output','-o', help = 'Output filename')
parser.add_argument('-z', action='store_true', default=False)
parser.add_argument('--extralabel', default='',help = 'Additional CMS label')

args = parser.parse_args()

with open(args.input) as jsonfile:
  js = json.load(jsonfile)

canv = ROOT.TCanvas(args.output,args.output)
pads = plot.OnePad()
pads[0].SetTicks(1,-1)
pads[0].SetLeftMargin(0.25)
pads[0].SetTicky(0)

#axis = ROOT.TH2F('axis', '',1,-2,2,7,0,7)
if args.z:
  xlo, xhi = 0, 2
else:
  xlo, xhi = -100, 300
axis = ROOT.TH2F('axis', '',1,xlo,xhi,10,0,10)
plot.Set(axis.GetYaxis(), LabelSize=0)
plot.Set(axis.GetXaxis(), Title = 'Best fit #mu')
axis.Draw()


cmb_band = ROOT.TBox()
plot.Set(cmb_band, FillColor=ROOT.kGreen)
if args.z:
  plot.DrawVerticalBand(pads[0],cmb_band,js['z']['val']-js['z']['ErrLo'],js['z']['val']+js['z']['ErrHi'])
else:
  plot.DrawVerticalBand(pads[0],cmb_band,js['r']['val']-js['r']['ErrLo'],js['r']['val']+js['r']['ErrHi'])

cmb_line = ROOT.TLine()
plot.Set(cmb_line,LineWidth=2)
if args.z:
  plot.DrawVerticalLine(pads[0],cmb_line,js['z']['val'])
else:
  plot.DrawVerticalLine(pads[0],cmb_line,js['r']['val'])

horizontal_line = ROOT.TLine()
plot.Set(horizontal_line,LineWidth=2,LineStyle=2)
plot.DrawHorizontalLine(pads[0],horizontal_line,3)

horizontal_line = ROOT.TLine()
plot.Set(horizontal_line,LineWidth=2,LineStyle=2)
plot.DrawHorizontalLine(pads[0],horizontal_line,6)


gr = ROOT.TGraphAsymmErrors(5)
plot.Set(gr, LineWidth=2, LineColor=ROOT.kRed)

#y_pos = 5.5
#x_text = -3.0
y_pos = 7.5
x_text = xlo - abs(xhi-xlo)*0.3
i=0
latex = ROOT.TLatex()
plot.Set(latex, TextAlign=12,TextSize=0.03)
latex.SetTextFont(42)
#order = ['r_ZH','r_WH','r_zerolep','r_onelep','r_twolep']
if args.z:
  order = ['z16', 'z17', 'z18', "z"]
  labels = []
else:
  order = ['r16', 'r17', 'r18', "r"]
labels = ["2016", "2017", "2018", "Run2"]

txt_dict = {}
for entry, label in zip(order, labels):
    txt_dict[entry] = '#splitline{%s}{#mu=%.2f#pm%.2f}'%(label,js[entry]['val'],(js[entry]['ErrHi']+js[entry]['ErrLo'])/2.)

for stre in order:
  gr.SetPoint(i,js[stre]['val'],y_pos)
  gr.SetPointError(i,js[stre]['ErrLo'],js[stre]['ErrHi'],0,0)
  latex.DrawLatex(x_text,y_pos,txt_dict[stre])

  i+=1
  y_pos -= 1.

gr.Draw("SAMEP")

pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

plot.DrawCMSLogo(pads[0],'CMS','%s'%args.extralabel,11,0.045,0.03,1.0,'',1.0)
#plot.DrawTitle(pads[0],'41.3 fb^{-1} (13 TeV)',3)
#plot.DrawTitle(pads[0],'39,6 fb^{-1} (13 TeV)',3)
plot.DrawTitle(pads[0],'137.6 fb^{-1} (13 TeV)',3)

latex.SetTextFont(42)
latex.SetTextSize(0.03)
#latex.DrawLatex(-0.82,6.1,"pp#rightarrow VH; H#rightarrow b#bar{b}")
#latex.DrawLatex(-0.82,5.7,"Combined #mu=%.1f#pm%.1f"%(js['z']['val'],js['z']['ErrHi']))
#latex.DrawLatex(-1.82,5.6,"pp#rightarrow VZ; H#rightarrow c#bar{c}")
if args.z:
  xpos = (xhi-xlo)/2 + (xhi-xlo) * 0.05
  latex.DrawLatex(xpos,9.4,"pp#rightarrow gg(Z#rightarrow c#bar{c})")
else:
  xpos = xlo + (xhi-xlo)/2 + (xhi-xlo) * 0.05
  latex.DrawLatex(xpos,9.4,"pp#rightarrow gg(H#rightarrow c#bar{c})")

#latex.DrawLatex(-1.82,5.2,"#mu=%.2f#pm%.2f(stat.)#pm%.2f(syst.)"%(js['z']['val'],(js['z']['StatHi']+js['z']['StatLo'])/2.,(js['z']['SystHi']+js['z']['SystLo'])/2.))
if args.z:
  latex.DrawLatex(xpos,8.8,"#mu=%.2f#pm%.2f(stat.+syst.)"%(js['z']['val'],(js['z']['ErrHi']+js['z']['ErrLo'])/2.))
else:
  latex.DrawLatex(xpos,8.8,"#mu=%.2f#pm%.2f(stat.+syst.)"%(js['r']['val'],(js['r']['ErrHi']+js['r']['ErrLo'])/2.))
#latex.DrawLatex(-1.8,5.2,"#mu=%.2f#pm%.2f(stat.+syst.)"%(js['z']['val'],(js['z']['ErrHi']+js['z']['ErrLo'])/2.))

canv.Print('.png')
canv.Print('.pdf')
