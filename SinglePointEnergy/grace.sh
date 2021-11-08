# make a panel plot
#arrange (3,3,.6,.6,.6,ON,ON,ON)
arrange (3,3,0.1,.5,.2,ON,OFF,OFF)
# chose the first panel
FOCUS G0
# I was hoping this line would allow me to change the axis limits, but it isn't working:
#each file has 6 columns
#s0 to s4
#READ NXY "file2.dat"
#s5 to s9
READ NXY "10_error.csv"
s0 linewidth 3 
s1 linewidth 3 
s2 linewidth 3 
s3 linewidth 3 
s4 linewidth 3 
s5 linewidth 3 
s0 line color 1
s0 linestyle 7
s1 line color 2
s2 linestyle 2
s2 line color 3
s3 linestyle 3
s3 line color 4
s3 linestyle 4
s4 line color 5
s4 linestyle 5
s5 line color 6
s5 linestyle 6
zeroyaxis on
zeroyaxis label "RCutCoulomb 10 \cE\C"
zeroyaxis label place opposite
s0 legend "Gross_DSP"
s1 legend "Gross_DSF"
s2 legend "Vlugt_DSP"
s3 legend "Vlugt_DSF"
s4 legend "Cassandra_DSP"
s5 legend "Cassandra_DSF"
xaxis tick place normal
xaxis TICKLABEL off
WORLD XMAX 0.3 
yaxis label "Relative Error"
# chose the first panel
FOCUS G1
READ NXY "11_error.csv"
s0 linewidth 3 
s1 linewidth 3 
s2 linewidth 3 
s3 linewidth 3 
s4 linewidth 3 
s5 linewidth 3 
s0 line color 1
s0 linestyle 7
s1 line color 2
s2 linestyle 2
s2 line color 3
s3 linestyle 3
s3 line color 4
s3 linestyle 4
s4 line color 5
s4 linestyle 5
s5 line color 6
s5 linestyle 6
zeroyaxis on
zeroyaxis label "RCutCoulomb 11 \cE\C"
zeroyaxis label place opposite
xaxis tick place normal
xaxis TICKLABEL off
WORLD XMAX 0.3 
yaxis label "Relative Error"
# chose the first panel
FOCUS G2
READ NXY "12_error.csv"
s0 linewidth 3 
s1 linewidth 3 
s2 linewidth 3 
s3 linewidth 3 
s4 linewidth 3 
s5 linewidth 3 
s0 line color 1
s0 linestyle 7
s1 line color 2
s2 linestyle 2
s2 line color 3
s3 linestyle 3
s3 line color 4
s3 linestyle 4
s4 line color 5
s4 linestyle 5
s5 line color 6
s5 linestyle 6
zeroyaxis on
zeroyaxis label "RCutCoulomb 12 \cE\C"
zeroyaxis label place opposite
xaxis label "Alpha"
xaxis tick place normal
WORLD XMAX 0.3 
yaxis label "Relative Error"
# chose the first panel
FOCUS G3
READ NXY "13_error.csv"
s0 linewidth 3 
s1 linewidth 3 
s2 linewidth 3 
s3 linewidth 3 
s4 linewidth 3 
s5 linewidth 3 
s0 line color 1
s0 linestyle 7
s1 line color 2
s2 linestyle 2
s2 line color 3
s3 linestyle 3
s3 line color 4
s3 linestyle 4
s4 line color 5
s4 linestyle 5
s5 line color 6
s5 linestyle 6
zeroyaxis on
zeroyaxis label "RCutCoulomb 13 \cE\C"
zeroyaxis label place opposite
xaxis TICKLABEL off
WORLD XMAX 0.3 
xaxis tick place normal
# chose the first panel
FOCUS G4
READ NXY "14_error.csv"
s0 linewidth 3 
s1 linewidth 3 
s2 linewidth 3 
s3 linewidth 3 
s4 linewidth 3 
s5 linewidth 3 
s0 line color 1
s0 linestyle 7
s1 line color 2
s2 linestyle 2
s2 line color 3
s3 linestyle 3
s3 line color 4
s3 linestyle 4
s4 line color 5
s4 linestyle 5
s5 line color 6
s5 linestyle 6
zeroyaxis on
zeroyaxis label "RCutCoulomb 14 \cE\C"
zeroyaxis label place opposite
xaxis TICKLABEL off
WORLD XMAX 0.3 
xaxis tick place normal
# chose the first panel
FOCUS G5
READ NXY "15_error.csv"
s0 linewidth 3 
s1 linewidth 3 
s2 linewidth 3 
s3 linewidth 3 
s4 linewidth 3 
s5 linewidth 3 
s0 line color 1
s0 linestyle 7
s1 line color 2
s2 linestyle 2
s2 line color 3
s3 linestyle 3
s3 line color 4
s3 linestyle 4
s4 line color 5
s4 linestyle 5
s5 line color 6
s5 linestyle 6
zeroyaxis on
zeroyaxis label "RCutCoulomb 15 \cE\C"
zeroyaxis label place opposite
xaxis label "Alpha"
WORLD XMAX 0.3 
xaxis tick place normal
KILL G6
KILL G7
KILL G8
PAGE SIZE 1200, 1000

