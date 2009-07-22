

import Gnuplot, Gnuplot.funcutils






#g.title('A simple example')


    
    
    



def PlotPopl(ratio,pop,p_size,data,pic):
    

    lines = open(data,"r").readlines()
    n_row_max = len(lines[0].split())
    xmin = float(lines[0].split()[0])
    xmax = float(lines[len(lines)-1].split()[0])
    xav = (xmin+xmax)*0.5



    g = Gnuplot.Gnuplot(debug=1)
    g('set term png ')
    g('set output "%s"'%pic)
    g('set xlabel "Time"')
    g('set size %f,%f'%(ratio[0],ratio[1]))
    g('set pointsize %f'%p_size)
    g('set multiplot')
    g('set yrange [0:1.01] ')
    g('set size 1.,1.')  
    
  
    n_row = 2

    for j in range(len(pop)):
        g(' set label "     %s" at %s,%f center font "Vera,30"  '%(pop[j],str(xav),0.8-0.1*j))

    for ny, nx  in [(ny,nx)
                            for ny in range(ratio[1]-1,-1,-1)
                            for nx in range(ratio[0])]:
                                
        if (n_row == 3):        g(' unset label ' )

        g('set title "000"')
        g('set origin %f,%f'%(nx,ny))
        g('plot "%s" using 1:%f title ""   '%(data,n_row))
        if (n_row == n_row_max):    break
        n_row += 1


    g('unset multiplot')




