import numpy as np
import matplotlib.pyplot as plt

def compute_slope(ax,xls,yls,x,y,slopes,R2,ipar,param_line,colors):
   #https://math.stackexchange.com/questions/3500898/understanding-the-least-squares-regression-formula
   #https://towardsdatascience.com/linear-regression-using-least-squares-a4c3456e8570
   # print('least-squares',len(xls),len(yls))

   # print(xls)
   # print(yls)
   logslope=True
   if logslope:
      X_mean1 = np.mean(xls)
      Y_mean1 = np.mean(yls)
      xls=np.log10(xls)
      yls=np.log10(yls)

   X_mean = np.mean(xls)
   Y_mean = np.mean(yls)
   X_min = np.min(xls)
   Y_min = np.min(yls)
   X_max = np.max(xls)
   Y_max = np.max(yls)

   num = 0
   den = 0
   for i in range(len(xls)):
      num += (xls[i] - X_mean)*(yls[i] - Y_mean)
      den += (xls[i] - X_mean)**2
   m = num / den
   c = Y_mean - m*X_mean

   Y_pred = m*xls + c

   rms=np.linalg.norm(yls-Y_pred, ord=2)

   rmsx=np.linalg.norm(xls-X_mean, ord=2)
   rmsy=np.linalg.norm(yls-Y_mean, ord=2)

   corr=num/(rmsx*rmsy)
   # print('corr',corr)

   if logslope:
      # ax.plot([10**(X_mean), 10**(max(xls))], [10**(X_mean*m+c), 10**(max(Y_pred))], color='black',alpha=0.5)
      line1,=ax.plot([10**(min(xls)), 10**(max(xls))], [10**(min(Y_pred)), 10**(max(Y_pred))],
                                 # color='black',
                                 color=colors[ipar+1],
                                 alpha=0.5,
                                 # label=str(m)
                                 )
      
      # param_line.append(line1)
      # ax.annotate('Slope '+"{:.1f}".format(m),xy=(X_mean1,Y_mean1))

      # print ('least-squares',x,y,m,c,10**(X_mean),10**(Y_mean))

      # ax.scatter(x=10**(X_mean),y=10**(Y_mean))
      # ax.scatter(x=10**(X_min),y=10**(Y_min))
      # ax.scatter(x=10**(X_max),y=10**(Y_max))

      # print(xls)
      # print(yls)


      # plt.text(0.8,0.9,
      # 'Slope={:.2g}\nR²{:.2g}'.format(float(m),float(rms)),
      # transform =ax.transAxes,fontsize=16,ha='center')
      slopes[ipar]=m
      R2[ipar]=corr#rms

      # ax.annotate('Slope={:.2g}\nR²={:.2g}'.format(float(m),float(rms)),
      # xy=(0.8,0.1),xycoords='axes fraction',ha='center')


      # plt.legend(('data', 'line-regression r={}'.format(r_value)), 'best')
      # legendslope=plt.legend(handles=[],labels=('Slope={:.2g}\nR²{:.2g}'.format(float(m),float(rms))), loc='best')
      # ax.add_artist(legendslope)

   else:   
      if (min(Y_pred)<0):
         ax.plot([X_mean, max(xls)], [X_mean*m+c, max(Y_pred)], color='black',
               #   alpha=0.5,
                  lw=lw) # predicted
      else:
         ax.plot([min(xls), max(xls)], [min(Y_pred), max(Y_pred)], color='black',
               #   alpha=0.5,
                  lw=lw) # predicted
   
      # ax.annotate("{:.1f}".format(m)+'*x'+"{:.3f}".format(c),xy=(X_mean,Y_mean))
      ax.annotate('Slope '+"{:.1f}".format(m),xy=(X_mean,Y_mean))

   print ('least-squares',x,y,m,c)#,10**(X_mean),10**(Y_mean))
   
   # print(X_mean,Y_mean,
   #       # len(xls),len(yls)
   #       ,min(xls),max(xls),min(Y_pred),max(Y_pred))

   print(min(Y_pred),max(Y_pred),Y_min,Y_max,10**(min(xls)),10**(max(xls)))   
   return(ax)


def plot_errors(domain_length,nx_list,l1_rel_error,l2_rel_error,linfty_rel_error):

    fig, ax = plt.subplots(layout="constrained")

    x = nx_list

    print(x)
    print(l1_rel_error)
    print(l2_rel_error)
    print(linfty_rel_error)


    plt.plot(x, l1_rel_error)
    plt.plot(x, l2_rel_error)
    plt.plot(x, linfty_rel_error)

    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)


    plt.savefig("errors.pdf")
    print("plot python")

    plt.close(fig)


# plot_bc2()
# print("test python pdi")
