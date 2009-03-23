

import Gnuplot, Gnuplot.funcutils






#g.title('A simple example')


    
    
    



def PlotPopl(ratio,pop,p_size,data,pic):
    
    f = open(data,"r")
    n_row_max = len(f.readline().split())
    f.close()
    
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

    for ny, nx  in [(ny,nx)
                            for ny in range(ratio[1]-1,-1,-1)
                            for nx in range(ratio[0])]:
                                
        if (n_row == n_row_max):    g(' set label "     %f" at 0.,0.5 center font "Vera,50"  '%pop)
        g('set title "000"')
        g('set origin %f,%f'%(nx,ny))
        g('plot "%s" using 1:%f title ""   '%(data,n_row))
        if (n_row == n_row_max):    break
        n_row += 1


    g('unset multiplot')




