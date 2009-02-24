from numpy import *

# If the package has been installed correctly, this should work:
import Gnuplot, Gnuplot.funcutils
import sys
import time




#g.title('A simple example')

def PlotPopl(population):
    
    
    
    g = Gnuplot.Gnuplot(debug=1)
    g('set term png ')
    g('set output "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/image.png"')
    g('set xlabel "Time"')
    g('set size 5.0,3.0')
    g('set pointsize 0.3')
    g('set multiplot')
    g('set yrange [0:1.01] ')
    g('set size 1.,1.')  
    
  


    
    g('set title "000"')
    g('set origin 0.,2')
    g('plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:2 title ""	')

    g('set title "00-1"')
    g('set origin 1,2')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:3 title ""	')

    g('set title "010"')
    g('set origin 2,2')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:4 title ""	')

    g('set title "100"')
    g('set origin 3,2')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:5 title ""')

    g('set title "001"')
    g('set origin 4,2')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:6 title ""')



    g('set title "00-2"')
    g('set origin 0.,1')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:7 title ""	')
    
    g('set title "01-1"')
    g('set origin 1,1')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:8 title ""	')
    
    g('set title "10-1"')
    g('set origin 2,1')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:9 title ""	')
    
    g('set title "020"')
    g('set origin 3,1')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:10 title ""	')
    
    g('set title "110"')
    g('set origin 4,1')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:11 title ""	')



    g('set title "200"')
    g('set origin 0.,0.')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:12 title ""	')

    g('set title "011"')
    g('set origin 1,0.')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:13 title ""	')

    g('set title "101"')
    g('set origin 2,0.')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:14 title ""	')

    g('set title "002"')
    g('set origin 3,0.')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:15 title ""	')
    

    
    g(' set label "     %f" at 0.,0.5 left font "Vera,50"  '%population)
    g('set title "Population Sum"')
    g('set origin 4,0.')
    g(' plot "/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt" using 1:16 title ""    ')





    g('unset multiplot')




