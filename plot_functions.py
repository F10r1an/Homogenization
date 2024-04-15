
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.stats as spstat 
import scipy.optimize as sciopt

# text rendering with LaTex
plt.rc('text', usetex=True)
plt.rc('font', family='times')

def plot_correction_tensor_regression_line(fig, pos, x, y, title, TOL=1.e-5, eps_l='', eps_s='', legend_loc=(-0.6, 1.0), col='red'):
    
    
    # use different bounds for optimizations if values are increasing or decreasing when converging to limit value
    # because log(x) has to be > 0
    

    def r_value(x0):
            lr = spstat.linregress(np.log10(x), np.log10(np.abs(np.array(y)-x0)))
            return lr.rvalue
    """
    # decreaseing case
    if y[-1] <= y[0]:  
        
        
        
        ropt = sciopt.minimize_scalar(r_value, method='bounded', bounds=[0.0, np.min(y) * (1. - TOL)] )
        y_opt = ropt.x
        #print(y_opt)
        lr = spstat.linregress(np.log10(x), np.log10(np.array(y)-y_opt))
        
        
    elif y[-1] > y[0]:  # increasing case
        def r_value(x0):
            lr = spstat.linregress(np.log10(x), np.log10(x0-np.array(y)))
            return lr.rvalue
        ropt = sciopt.minimize_scalar(r_value, method='bounded', bounds=[np.max(y)*(1+TOL),2*np.max(y) ] )
        y_opt = ropt.x
        lr = spstat.linregress(np.log10(x), np.log10(y_opt-np.array(y)))
    """
    ropt = sciopt.minimize_scalar(r_value, method='bounded', bounds=[0.0,2*np.max(y) ] )
    y_opt = ropt.x
    lr = spstat.linregress(np.log10(x), np.log10(np.abs(y_opt-np.array(y))))
        
    
    r_value = lr.rvalue
    p_value = lr.pvalue
    slope = lr.slope
    intercept = lr.intercept
    stderr = lr.stderr
        
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.set_xlabel(r'$log_{10}$ (\# cells)', fontsize=8)
    if title == 'Diffusion':
        sub.set_ylabel(r'$log_{10}$ ($D_{11} / \theta / D_0$) ', fontsize=8)
    elif title == 'Permittivity':
        sub.set_ylabel(r'$log_{10}$ ($(E_{11}-\varepsilon_s) /( \varepsilon_l-\varepsilon_s)$)', fontsize=8)
    elif title == 'Mobility':
        sub.set_ylabel(r'$log_{10}$ ($M_{11}/ \theta / D_0$)', fontsize=8)
    
    #sub.text(0., 1.05, r'$\varepsilon_l$/$\varepsilon_s=$'+str(eps_l), transform=sub.transAxes, fontsize=10)
    #sub.set_title(title, pad=-50)
    
    if y[-1] < y[0]:
        if title == 'Diffusion':
            label =  ( r'$r^2$ value: ' + str(np.round(r_value*r_value,7))
                      + '\np value: ' + str(np.round(p_value,7))
                      + '\nstderr: ' + str(np.round(stderr, 7))
                      +'\nslope: ' + str(np.round(lr.slope, 4))
                      +'\nintercept: ' + str(np.round(lr.intercept, 4))
                      +'\noffset: ' + str(np.round(y_opt,4))
                       )
        else: 
        #
            label = (#'$ \epsilon_l=$'+str(eps_l)+'\n'
                      #+'\n 
                      '$\epsilon_s = $'+str(eps_s)+'\n'
                      +'$r^2$ value: ' + str(np.round(r_value*r_value,7))
                      + '\np value: ' + str(np.round(p_value,7))
                      + '\nstderr: ' + str(np.round(stderr, 7))
                      
                      #+'\nMin. value: ' + str(np.round(np.min(y), 4))
                      +'\nslope: ' + str(np.round(lr.slope, 4))
                      +'\nintercept: ' + str(np.round(lr.intercept, 4)) 
                      +'\noffset: ' + str(np.round(y_opt,4))
                      )
        sub.plot(np.log10(x), np.log10(np.abs(np.array(y)-y_opt)), "+", color=col, zorder=10)
        xfit = np.linspace(np.log10(np.min(x)), np.log10(np.max(x)), 1000)
        #yfit = -lr.intercept * np.power(xfit, lr.slope) 
        yfit = lr.slope * xfit + lr.intercept
        sub.plot(xfit, yfit, '--', color='black', 
                      label= label
                      )
    
    elif y[-1] > y[0]: 
        sub.plot(np.log10(x), np.log10(np.abs(y_opt-np.array(y))), "+", color=col, zorder=10)
        xfit = np.linspace(np.log10(np.min(x)), np.log10(np.max(x)), 1000)
        #yfit = -lr.intercept * np.power(xfit, lr.slope) 
        yfit = lr.slope * xfit + lr.intercept
        label =  (
        '$ \epsilon_l = $'+str(eps_l)+'\n'

        + '$\epsilon_s = $'+str(eps_s)+'\n'
                      +'$r^2$ value: ' + str(np.round(r_value*r_value,7))
                      + '\np value: ' + str(np.round(p_value,7))
                      + '\nstderr: ' + str(np.round(stderr, 7))
                      +'\nslope: ' + str(np.round(lr.slope, 4))
                      +'\nintercept: ' + str(np.round(lr.intercept, 4))
                      +'\noffset: ' + str(np.round(y_opt,4))
                       )
        sub.plot(xfit, yfit, '--', color='black', 
                      label=label
                      )
    
    sub.tick_params('both', labelsize=8)
    
    leg = sub.legend(loc=legend_loc, fontsize=8, facecolor=None, shadow=None, )
    leg.get_frame().set_linewidth(0.0)
    leg.set_zorder(-100)
    
    #sub.set_xscale('log')
    #sub.set_yscale('log')
    
    return y_opt

def plot_supp1b(fig, pos, x, y):
    """
    plot similar to plot_single_correction_tensor_vs_number_of_cells() but
    with broken y-axis
    """
    
    x0 = pos[0]
    dx = pos[2]
   
   
    
    d1 = 0.1
    d2 = 0.05  # whitespace
    d3 = 0.8
    y1 = pos[1]
    y2 = y1 + (d1+d2) * pos[3]
    dy1 = pos[3] * d1
    dy2 = pos[3] * d3
    
    
    # upper plot contains data points
    ax_upper = fig.add_axes([x0, y2, dx, dy2])
    ax_upper.set_xlim([-5000, 90000])
    ax_upper.set_ylim([0.362, 0.386])
    
    ax_upper.plot(x, y, '+', label='numerical results\n(this paper)', color='green')
    ax_upper.plot([144], [0.3833], 'x', label='numerical result \n (Auriault \& \n Lewandowska 1996)', color='black')
    ax_upper.plot([-5000, 90000], [0.3671, 0.3671], '--', label='experiement\n (Auriault \& \n Lewandowska 1996)', color='black')
    
    ax_upper.plot([-5000, 90000], [0.3641, 0.3641], '--', label='limit', color='darkorchid')
    
    ax_upper.spines['top'].set_visible(False)
    ax_upper.spines['right'].set_visible(False)
    ax_upper.spines['bottom'].set_visible(False)
    
    ax_upper.tick_params(labelsize=8)   
    
    
    leg=ax_upper.legend(fontsize=7, loc=(0.25,0.3), frameon=False)
    ax_upper.set_xticks([])
    ax_upper.set_ylabel(r'$D_{11} / \theta / D_0$ ',fontsize=8)  
    
    
    
    ############################
    ax_lower = fig.add_axes([x0, y1, dx, dy1])

    ax_lower.set_xlim((ax_upper.get_xlim()))
    ax_lower.set_facecolor('none')
    ax_lower.spines['top'].set_visible(False)
    ax_lower.spines['right'].set_visible(False)
    ax_lower.set_yticks([0])

    ax_lower.set_xticks(ticks=[0, 25000, 50000, 75000])
    ax_lower.set_xlabel('\# cells',fontsize=8)
    ax_lower.tick_params(labelsize=8)   
    
    #######################
    xl = 3000
    x0 = pos[0] - pos[2] / 95000 * xl  # shift to left
    dx = pos[2] * (95000+xl)/(95000)  # compensate width
    
    ax_label = fig.add_axes([x0, y1, dx, dy1*2.5])
    ax_label.set_ylim([0, 0.25])
    ax_label.set_xlim([-5000 -xl, 90000])
    
    ax_label.tick_params(labelsize=8)   
    
    ax_label.set_facecolor('none')
    
    # broken axis
    dy = 0.02
    ax_label.plot([-5000 - xl, -5000 + xl],[0.1+dy, 0.1-dy], c='k', lw=1.)
    ax_label.plot([-5000 - xl, -5000 + xl],[0.15+dy, 0.15-dy], c='k', lw=1.)
    # measure offset
    #ax_label.plot([80000,80000],[0, 0.23], c='k' , lw=1.)
    
    ax_label.set_yticks([])
    ax_label.set_xticks(ticks=[0, 25000, 50000, 75000])   
    ax_label.set_xticks([]) 
    ax_label.spines['left'].set_visible(False)
    ax_label.spines['right'].set_visible(False)
    ax_label.spines['bottom'].set_visible(False)
    ax_label.spines['top'].set_visible(False)
    
def plot_single_correction_tensor_vs_number_of_cells(fig, pos, x, y, title='', options=['plot_diagonals'], component='11'):
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.set_xlabel('\# cells',fontsize=8)
    if component == '11':
        sub.set_ylabel(r'$D_{11} / \theta / D_0$ ~[a. u.]',fontsize=8)  
    
    sub.plot(x, y, '+', label='Numerical results\n(this paper)', color='green')
    
    if 'compare_auriault' in options:  
    
        sub.plot([144], [0.3833], 'x', label='Numerical result \n (Auriault \& \n Lewandowska 1996)', color='black')
        sub.plot([144, np.max(x)], [0.3671, 0.3671], '--', label=' Experiement\n (Auriault \& \n Lewandowska 1996)', color='black')
    
    sub.tick_params(labelsize=8)
    
    leg=sub.legend(fontsize=8, loc=(0.25,0.35))
    leg.get_frame().set_linewidth(0.0)
       
def plot_correction_tensors_vs_number_of_cells(fig, pos1, pos2, x1, y1, x2, y2, title='', options=['plot_diagonals']):
    sub1 = fig.add_axes(pos1)
    sub1.spines['top'].set_visible(False)
    sub1.spines['right'].set_visible(False)
    sub1.set_xlabel('\# cells')
    sub1.set_ylabel(r'$D_{ij} / \theta / D_0$')
    
    if 'plot_all_components' in options:
        sub1.plot(x1, y1[0], '+', label='$D_{11}$')
        sub1.plot(x1, y1[1], 'x', label='$D_{12}$')
        sub1.plot(x1, y1[2], '+', label='$D_{21}$')
        sub1.plot(x1, y1[3], 'x', label='$D_{22}$')
    if 'plot_diagonals' in options:   
        sub1.plot(x1, np.array(y1[0]), '+', label='$D_{11}$')
        sub1.plot(x1, np.array(y1[1]), 'x', label='$D_{22}$')
    if 'plot_only_11_diffusion' in options:
        sub1.plot(x1, y1[0], '+', label='$D_{11}$')
    if 'compare_auriault' in options:
        sub1.plot([144], [0.3833], 'x', label='Auriault \& \n Lewandowska', color='black')
        sub1.plot([144, np.max(x1)], [0.3671, 0.3671], '--', label='Experiement', color='black')
    if 'plot_straight_channel_theory' in options:
        sub1.plot([np.min(x1), np.max(x1)], [1., 1.], '--', label='theory', color='black')
    if 'plot_diagonal_channel_theory' in options:
        sub1.plot([np.min(x1), np.max(x1)], [.5, .5], '--', label='theory', color='black')
        

    
    #sub1.plot([20,200], [0.3833, 0.3833], '--', color='black', label='Auriault \&\n Lewandowska 1996\n (144 cells)')
    #sub1.plot(res_list, d12, 'x', label='$D_{12}$')
    
    #sub1.plot(res_list, d21, 'x', label='$D_{21}$')

    #sub1.set_xticklabels([str(np.square(res_list[0])*2), '', '', str(np.square(res_list[3])*2),
    #                      '', '', str(np.square(res_list[6])*2), '', '', str(np.square(res_list[9])*2)])
    #sub1.set_xticklabels([str(np.square(res_list[0])*2), '', str(np.square(res_list[2])*2), '', str(np.square(res_list[4])*2)])
    
    #sub1.set_xscale('log')
    #sub1.set_yscale('log')
    sub1.legend()
    sub1.set_title('Diffusion')

    sub2 = fig.add_axes(pos2)
    sub2.spines['top'].set_visible(False)
    sub2.spines['right'].set_visible(False)
    sub2.set_xlabel('\# cells')
    sub2.set_ylabel('$E_{ij} ~[a. u.]$')
    
    if 'plot_all_components' in options:
        sub2.plot(x2, y2[0], '+', label='$E_{11}$')
        sub2.plot(x2, y2[1], 'x', label='$E_{12}$')
        sub2.plot(x2, y2[2], '+', label='$E_{21}$')
        sub2.plot(x2, y2[3], 'x', label='$E_{22}$')
    
    if 'plot_diagonals' in options or 'plot_permittivity_diagonals' in options:   
        sub2.plot(x2, y2[0], '+', label='$E_{11}$')
        sub2.plot(x2, y2[1], 'x', label='$E_{22}$')
    if 'plot_straight_channel_theory' in options:
        sub2.plot([np.min(x2), np.max(x2)], [1.333, 1.333], '--', label='theory', color='black')
        sub2.plot([np.min(x2), np.max(x2)], [1.5, 1.5], '--', color='black')
    if 'plot_diagonal_channel_theory' in options:
        sub2.plot([np.min(x2), np.max(x2)], [1.4167, 1.4167], '--', label='theory', color='black')
        sub2.plot([np.min(x2), np.max(x2)], [0.0833, 0.0833], '--', color='black')
        
    #sub2.plot([20, 200], [1.5, 1.5], '--', color='black', label='theory')
    #sub2.plot([20, 200], [1.333, 1.333], '--', color='black')
    #sub2.plot(res_list, p12, 'x', label='$E_{12}$')
    
    #sub2.plot(res_list, p21, 'x', label='$E_{21}$')

    #sub2.set_xticks(res_list)
    #sub2.set_xticklabels([str(np.square(res_list[0])*2), '', '', str(np.square(res_list[3])*2),
    #                      '', '', str(np.square(res_list[6])*2), '', '', str(np.square(res_list[9])*2)])pl
    #sub2.set_xticklabels([str(np.square(res_list[0])*2), '', str(np.square(res_list[2])*2), '', str(np.square(res_list[4])*2)])
    
    #sub2.set_xscale('log')
    #sub2.set_yscale('log')
    sub2.legend()
    sub2.set_title('Permittivity')

    if not title == '':
        fig.suptitle(title)
        
def plot_box_cell_is_dimensionless(fig, pos1, pos2, x1, y1, x2, y2, title=''):
    sub1 = fig.add_axes(pos1)
    sub1.spines['top'].set_visible(False)
    sub1.spines['right'].set_visible(False)
    sub1.set_xlabel('\# cells')
    sub1.set_ylabel(r'$D_{11}/D_0$')
    

    sub1.plot(x1[0], y1[0], '+', label='one box')
    sub1.plot(x1[1], y1[1], 'x', label='four boxes')
    sub1.plot(x1[2], y1[2], '1', label='crossing channels')
    sub1.plot(x1[3], y1[3], '2', label='diagonal box')
    #sub1.set_xscale('log')


    #sub1.legend()
    sub1.set_title('Diffusion')

    sub2 = fig.add_axes(pos2)
    sub2.spines['top'].set_visible(False)
    sub2.spines['right'].set_visible(False)
    sub2.set_xlabel('\# cells')
    sub2.set_ylabel('$E_{11}$')
    #sub2.set_xscale('log')
     
    sub2.plot(x2[0], y2[0], '+', label='one box')
    sub2.plot(x2[1], y2[1], 'x', label='four boxes')
    sub2.plot(x2[2], y2[2], '1', label='crossing channels')
    sub2.plot(x2[3], y2[3], '2', label='diagonal box')

    sub2.legend(loc=(1.0, .4))
    sub2.set_title('Permittivity')
    
    sub2.set_xscale('log')
    sub2.set_yscale('log')

    if not title == '':
        fig.suptitle(title)
        
def plot_scale_independence_permittivity(fig, pos, x, y, component='11', title=''):
    sub1 = fig.add_axes(pos)
    sub1.spines['top'].set_visible(False)
    sub1.spines['right'].set_visible(False)
    sub1.set_xlabel(r'$\varepsilon_s / \varepsilon_l$')
    # sub1.set_ylabel(r'($E_{11} - \varepsilon_l$) / ($\varepsilon_l - \varepsilon_s$)')
    if component == '11':
        sub1.set_ylabel(r'$E_{11} / \varepsilon_l$')
    elif component == '22':
        sub1.set_ylabel(r'$E_{22} / \varepsilon_l$')
   
    sub1.plot(x[0], y[0], 'x', label=r'$\varepsilon_l = 1$')
    sub1.plot(x[1], y[1], '+', label=r'$\varepsilon_l = 3$')
    sub1.plot(x[2], y[2], '1', label=r'$\varepsilon_l = 9$')
    sub1.plot(x[3], y[3], '2', label=r'$\varepsilon_l = 27$')
    sub1.plot(x[4], y[4], '3', label=r'$\varepsilon_l = 81$')
    
    sub1.set_xscale('log')
 
    
    sub1.legend()
    if title=='':
        sub1.set_title(r'Permittivity tensor scales dimensionless')
    else:
        sub1.set_title(title)
    
def plot_pourosity_vs_correction_tensors(fig, pos1, pos2, x1, y1, x2, y2, colors, title='', leg_loc=None, make_ylabel=True):
    #sub1 = fig.add_axes(pos1)
    sub1 = fig.add_axes(pos2)
    sub1.spines['top'].set_visible(False)
    sub1.spines['right'].set_visible(False)
    

    

    sub1.set_xlim([-0.05, 1.05])   
    sub1.set_ylim([-0.05, 1.05])   
    
    sub1.plot(x1[0], y1[0], '+', label='$D_{11}$', color='grey')
    sub1.plot([0.,1.],[0.,1.], '--', color='black') 

    sub1.set_xticks([0.,0.5,1.])
    # sub1.set_xticklabels([0,'',1])
    sub1.set_xticklabels(['','',''])
    
    sub1.set_yticks([0., 0.5, 1.])
    if make_ylabel == True:
        sub1.set_ylabel(r'$D_{ii} / D_0$', labelpad=0, fontsize=8)        
        sub1.set_yticklabels([0, '', 1])
    else: 
        sub1.set_ylabel('')        
        sub1.set_yticklabels(['', '', ''])
    #sub1.set_xlabel(r'Porosity $\theta$', fontsize=8, labelpad=0)
    sub1.tick_params(labelsize=8)
    
    #sub2 = fig.add_axes(pos2)
    sub2 = fig.add_axes(pos1)
    sub2.spines['top'].set_visible(False)
    sub2.spines['right'].set_visible(False)
    sub2.tick_params(labelsize=8)


    sub2.set_xlim([-0.05, 1.05])  
    sub2.set_ylim([-0.05, 1.05])  
    sub2.set_xticks([0.,0.5,1.])
    sub2.set_yticks([0., 0.5, 1.])
    #sub2.set_xticklabels(['','',''])
    sub2.set_xticklabels([0,'',1])
    sub2.set_xlabel('')
    sub2.set_xlabel(r'Porosity $\theta$', fontsize=8, labelpad=0)
    
    sub2.plot(x2[0], y2[0], 'x', label=r'$\varepsilon_l / \varepsilon_s=100.$', c=colors[0])
    sub2.plot(x2[1], y2[1], 'x', label=r'$\varepsilon_l / \varepsilon_s=10.$', c=colors[1])
    sub2.plot(x2[2], y2[2], 'x', label=r'$\varepsilon_l / \varepsilon_s=2.$', c=colors[2])
    sub2.plot(x2[3], y2[3], 'x', label=r'$\varepsilon_l / \varepsilon_s=0.5$', c=colors[3])
    sub2.plot(x2[4], y2[4], 'x', label=r'$\varepsilon_l / \varepsilon_s=0.1$', c=colors[4])
    sub2.plot(x2[5], y2[5], 'x', label=r'$\varepsilon_l / \varepsilon_s=0.01$', c=colors[5])
    
    sub2.plot([0.,1.],[0.,1.], '--', color='black')
    
    
    
    if make_ylabel == True:
        sub2.set_ylabel(r'($E_{ii}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$)', fontsize=8, labelpad=0)
        sub2.set_yticklabels([0, '', 1])
    else: 
        sub2.set_ylabel('')        
        sub2.set_yticklabels(['', '', ''])

    #sub2.set_xscale('log')
    if leg_loc == None:
        pass
    else:
        sub2.legend(loc=leg_loc, ncol=6, fontsize=8, frameon=False)
    
    #sub2.set_title(title)

    # fig.suptitle(title)
    
def plot_pourosity_vs_correction_tensors_with_mobility(fig, pos1, pos2, pos3, x1, y1, x2, y2, x3, y3, title='', pos_legend=None, component='11'):
    
    ################ DIFFUSION
    sub1 = fig.add_axes(pos1)
    sub1.spines['top'].set_visible(False)
    sub1.spines['right'].set_visible(False)
    sub1.set_xlabel(r'Porosity $\theta$')
    if component == '11':
        sub1.set_ylabel(r'$D_{11} / D_0$ ~[a. u.]')
    elif component == '22':
        sub1.set_ylabel(r'$D_{22} / D_0$ ~[a. u.]')
    elif component == '21':
        sub1.set_ylabel(r'$D_{21} / D_0$ ~[a. u.]')
    elif component == '12':
        sub1.set_ylabel(r'$D_{12} / D_0$ ~[a. u.]')
    sub1.set_xlim([-0.05, 1.05])   
    sub1.set_ylim([-0.05, 1.05])   
    
    sub1.plot(x1[0], y1[0], '+', label='$D_{11}$', color='grey')
    sub1.plot([0.,1.],[0.,1.], '--', color='black') 

    #sub1.set_xscale('log')
    #sub1.legend()
    sub1.set_title('Diffusion')
    
    #################### PERMITTIVITY

    sub2 = fig.add_axes(pos2)
    sub2.spines['top'].set_visible(False)
    sub2.spines['right'].set_visible(False)
    sub2.set_xlabel(r'Porosity $\theta$')
    if component == '11':
        sub2.set_ylabel(r'($E_{11}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$) [a. u.]')
    elif component == '22':
        sub2.set_ylabel(r'($E_{22}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$) [a. u.]')
    elif component == '12':
        sub2.set_ylabel(r'($E_{12}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$) [a. u.]')
    elif component == '21':
        sub2.set_ylabel(r'($E_{21}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$) [a. u.]')
    sub2.set_xlim([-0.05, 1.05])  
    sub2.set_ylim([-0.05, 1.05])  
   
    
    
    sub2.plot(x2[0], y2[0], 'x', label=r'$\varepsilon_l / \varepsilon_s=100.$')
    sub2.plot(x2[1], y2[1], 'x', label=r'$\varepsilon_l / \varepsilon_s=10.$')
    sub2.plot(x2[2], y2[2], '+', label=r'$\varepsilon_l / \varepsilon_s=2.$')
    sub2.plot(x2[3], y2[3], 'x', label=r'$\varepsilon_l / \varepsilon_s=0.5$')
    sub2.plot(x2[4], y2[4], 'x', label=r'$\varepsilon_l / \varepsilon_s=0.1$')
    sub2.plot(x2[5], y2[5], 'x', label=r'$\varepsilon_l / \varepsilon_s=0.01$')
    
    sub2.plot([0.,1.],[0.,1.], '--', color='black')

    #sub2.set_xscale('log')
    if pos_legend == None:
        sub2.legend(loc=(1.1, 0.2))
    else: 
        sub2.legend(loc=pos_legend)
        
    sub2.set_title(r'Permittivity')
    
    
    ##################### MOBILITY
    sub3 = fig.add_axes(pos3)
    sub3.spines['top'].set_visible(False)
    sub3.spines['right'].set_visible(False)
    sub3.set_xlabel(r'Porosity $\theta$')
    if component == '11':
        sub3.set_ylabel(r'$M_{11} / D_0$ ~[a. u.]')
    elif component == '22':
        sub3.set_ylabel(r'$M_{22} / D_0$ ~[a. u.]')
    elif component == '12':
        sub3.set_ylabel(r'$M_{12} / D_0$ ~[a. u.]')
    elif component == '21':
        sub3.set_ylabel(r'$M_{21} / D_0$ ~[a. u.]')
    sub3.set_xlim([-0.05, 1.05])  
    sub3.set_ylim([-0.05, 1.05])  

    fig.suptitle(title)
    
    sub3.plot(x3[0], y3[0], 'x')
    sub3.plot(x3[1], y3[1], 'x')
    sub3.plot(x3[2], y3[2], '+')
    sub3.plot(x3[3], y3[3], 'x')
    sub3.plot(x3[4], y3[4], 'x')
    sub3.plot(x3[5], y3[5], 'x')
    
    sub3.plot([0.,1.],[0.,1.], '--', color='black')
    sub3.set_title(r'Mobility')
    
    if title == 'Straight channel 2D':
        sub1.plot(straight_channel_theory_diffusion('11')[0], straight_channel_theory_diffusion('11')[1], color='lightgrey', linestyle='-', zorder=-100, lw=5.)
        sub2.plot(straight_channel_theory_permittivity(100., 1., component='11')[0],
                  straight_channel_theory_permittivity(100., 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100, lw=5.)
        sub2.plot(straight_channel_theory_permittivity(10., 1., component='11')[0],
                  straight_channel_theory_permittivity(10., 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(2., 1., component='11')[0],
                  straight_channel_theory_permittivity(2., 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(0.5, 1., component='11')[0],
                  straight_channel_theory_permittivity(0.5, 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(0.1, 1., component='11')[0],
                  straight_channel_theory_permittivity(0.1, 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(0.01, 1., component='11')[0],
                  straight_channel_theory_permittivity(0.01, 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
    #straight_channel_theory_permittivity(component='')
    #diagonal_channel_theory_permittivity(eps_l, eps_s, component='')
    #diagonal_channel_theory_diffusion(component='')
    if title == 'Straight channel 2D comp22':
        sub1.plot(straight_channel_theory_diffusion('22')[0], straight_channel_theory_diffusion('22')[1], color='lightgrey', linestyle='--', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(100., 1., component='22')[0],
                  straight_channel_theory_permittivity(100., 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(10., 1., component='22')[0],
                  straight_channel_theory_permittivity(10., 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(2., 1., component='22')[0],
                  straight_channel_theory_permittivity(2., 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(0.5, 1., component='22')[0],
                  straight_channel_theory_permittivity(0.5, 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(0.1, 1., component='22')[0],
                  straight_channel_theory_permittivity(0.1, 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(straight_channel_theory_permittivity(0.01, 1., component='22')[0],
                  straight_channel_theory_permittivity(0.01, 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
    if title == 'Diagonal channel 2D':
        sub1.plot(diagonal_channel_theory_diffusion('11')[0], diagonal_channel_theory_diffusion('11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(100., 1., component='11')[0],
                  diagonal_channel_theory_permittivity(100., 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(10., 1., component='11')[0],
                  diagonal_channel_theory_permittivity(10., 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(2., 1., component='11')[0],
                  diagonal_channel_theory_permittivity(2., 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.5, 1., component='11')[0],
                  diagonal_channel_theory_permittivity(0.5, 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.1, 1., component='11')[0],
                  diagonal_channel_theory_permittivity(0.1, 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.01, 1., component='11')[0],
                  diagonal_channel_theory_permittivity(0.01, 1., component='11')[1], color='lightgrey', linestyle='-', zorder=-100)
    if title == 'Diagonal channel 2D comp22':
        sub1.plot(diagonal_channel_theory_diffusion('22')[0], diagonal_channel_theory_diffusion('22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(100., 1., component='22')[0],
                  diagonal_channel_theory_permittivity(100., 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(10., 1., component='22')[0],
                  diagonal_channel_theory_permittivity(10., 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(2., 1., component='22')[0],
                  diagonal_channel_theory_permittivity(2., 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.5, 1., component='22')[0],
                  diagonal_channel_theory_permittivity(0.5, 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.1, 1., component='22')[0],
                  diagonal_channel_theory_permittivity(0.1, 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.01, 1., component='22')[0],
                  diagonal_channel_theory_permittivity(0.01, 1., component='22')[1], color='lightgrey', linestyle='-', zorder=-100)
    
    if title == 'Diagonal channel 2D comp12':
        sub1.plot(diagonal_channel_theory_diffusion('12')[0], diagonal_channel_theory_diffusion('12')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(100., 1., component='12')[0],
                  diagonal_channel_theory_permittivity(100., 1., component='12')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(10., 1., component='12')[0],
                  diagonal_channel_theory_permittivity(10., 1., component='12')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(2., 1., component='12')[0],
                  diagonal_channel_theory_permittivity(2., 1., component='12')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.5, 1., component='12')[0],
                  diagonal_channel_theory_permittivity(0.5, 1., component='12')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.1, 1., component='12')[0],
                  diagonal_channel_theory_permittivity(0.1, 1., component='12')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.01, 1., component='12')[0],
                  diagonal_channel_theory_permittivity(0.01, 1., component='12')[1], color='lightgrey', linestyle='-', zorder=-100)
        
        sub2.set_ylim(-1.1, 2.1)
        sub3.set_ylim(-1.1, 2.1)
    
    if title == 'Diagonal channel 2D comp21':
        sub1.plot(diagonal_channel_theory_diffusion('21')[0], diagonal_channel_theory_diffusion('21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(100., 1., component='21')[0],
                  diagonal_channel_theory_permittivity(100., 1., component='21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(10., 1., component='21')[0],
                  diagonal_channel_theory_permittivity(10., 1., component='21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(2., 1., component='21')[0],
                  diagonal_channel_theory_permittivity(2., 1., component='21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.5, 1., component='21')[0],
                  diagonal_channel_theory_permittivity(0.5, 1., component='21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.1, 1., component='21')[0],
                  diagonal_channel_theory_permittivity(0.1, 1., component='21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.plot(diagonal_channel_theory_permittivity(0.01, 1., component='21')[0],
                  diagonal_channel_theory_permittivity(0.01, 1., component='21')[1], color='lightgrey', linestyle='-', zorder=-100)
        sub2.set_ylim(-1.1, 2.1)
        sub3.set_ylim(-1.1, 2.1)       
        
def plot_different_box_reference_cells(fig, pos):
    #height = fig.get_figheight()
    #width = fig.get_figwidth()
    scale=18.
    bs = 8.
    offset = 0.5*scale-0.5*bs
    n=5
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.spines['bottom'].set_visible(False)
    sub.spines['left'].set_visible(False)
    sub.set_xlim(0,n*scale)
    sub.set_ylim(0,n*(scale))
    sub.set_xticks([])
    sub.set_yticks([])
    sub.fill_between([0,n*scale],[0,0],[n*(scale),n*(scale)],color='steelblue', label='Liquid phase')
    
    for x in range(n):
        for y in range(n):
            if x == 1 and y==1:
                sub.fill_between([x*scale+offset, x*scale+offset+bs], 
                [y*scale+offset, y*scale+offset], [y*scale+offset+bs, y*scale+offset+bs], 
                color='darkgrey', label='Solid phase')
            else:
                sub.fill_between([x*scale+offset, x*scale+offset+bs], 
                [y*scale+offset, y*scale+offset], [y*scale+offset+bs, y*scale+offset+bs], 
                color='darkgrey')
    
    # ref cell 1
    c1 = 'sandybrown'
    c2 = 'darkorchid'
    c3 = 'turquoise'
    c4 = 'limegreen'
    
    sub.plot([2.5*scale, 3.5*scale], [0.5*scale, 0.5*scale], color=c2, label='Crossing channels 2D' )
    sub.plot([2.5*scale, 3.5*scale], [1.5*scale, 1.5*scale], color=c2 )
    sub.plot([2.5*scale, 2.5*scale], [0.5*scale, 1.5*scale],  color=c2 )
    sub.plot([3.5*scale, 3.5*scale], [0.5*scale, 1.5*scale], color=c2 )
    
    sub.plot([3.5*scale, 4.5*scale], [2.5*scale, 3.5*scale], color=c3, label='Diagonal box in a box 2D' )
    sub.plot([4.5*scale, 3.5*scale], [3.5*scale, 4.5*scale], color=c3 )
    sub.plot([3.5*scale, 2.5*scale], [4.5*scale, 3.5*scale], color=c3 )
    sub.plot([2.5*scale, 3.5*scale], [3.5*scale, 2.5*scale], color=c3 )
    
    sub.plot([1*scale, 3*scale], [1*scale, 1*scale], color=c4, label='Four boxes in a box 2D' )
    sub.plot([1*scale, 3*scale], [3*scale, 3*scale], color=c4 )
    sub.plot([1*scale, 1*scale], [1*scale, 3*scale], color=c4 )
    sub.plot([3*scale, 3*scale], [1*scale, 3*scale], color=c4 )
    
    sub.plot([1*scale, 2*scale], [1*scale, 1*scale], '--', color=c1 )
    sub.plot([1*scale, 2*scale], [2*scale, 2*scale], color=c1, label='Box in a box 2D' )
    sub.plot([1*scale, 1*scale], [1*scale, 2*scale], '--', color=c1 )
    sub.plot([2*scale, 2*scale], [1*scale, 2*scale], color=c1 )
    
    #leg = plt.legend(loc=(1.1, 0.4))
    #leg.get_frame().set_linewidth(0.0)
    
def plot_different_channel_reference_cells(fig, pos):
    scale=9.
    cs = 4.5
    offset = 0.5*scale-0.5*cs
    n=3
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.spines['bottom'].set_visible(False)
    sub.spines['left'].set_visible(False)
    sub.set_xlim(0,n*scale)
    sub.set_ylim(0,n*(scale))
    sub.set_xticks([])
    sub.set_yticks([])
    sub.fill_between([0,n*scale],[0,0],[n*(scale),n*(scale)],color='steelblue', label='Liquid phase')
    
    
    for y in range(n):
        if y==1:
            sub.fill_between([0, n*scale], 
            [y*scale+offset, y*scale+offset], [y*scale+offset+cs, y*scale+offset+cs], 
            color='darkgrey', label='Solid phase')
        else:
            sub.fill_between([0, n*scale], 
            [y*scale+offset, y*scale+offset], [y*scale+offset+cs, y*scale+offset+cs], 
            color='darkgrey')

    # ref cell 1
    c1 = 'crimson'
    c2 = 'gold'
    
    
    sub.plot([0.2*scale, 1.2*scale], [0.5*scale-1.5/9*scale, 0.5*scale-1.5/9*scale], color=c1, label='Straight channel 2D' )
    sub.plot([0.2*scale, 1.2*scale], [1.5*scale-1.5/9*scale, 1.5*scale-1.5/9*scale], color=c1 )
    sub.plot([0.2*scale, 0.2*scale], [0.5*scale-1.5/9*scale, 1.5*scale-1.5/9*scale],  color=c1 )
    sub.plot([1.2*scale, 1.2*scale], [0.5*scale-1.5/9*scale, 1.5*scale-1.5/9*scale], color=c1 )
    
    sub.plot([1.8*scale, 2.8*scale], [.5*scale, 1.5*scale], color=c2, label='Diagonal channel 2D' )
    sub.plot([2.8*scale, 1.8*scale], [1.5*scale, 2.5*scale], color=c2 )
    sub.plot([1.8*scale, 0.8*scale], [2.5*scale, 1.5*scale], color=c2 )
    sub.plot([0.8*scale, 1.8*scale], [1.5*scale, 0.5*scale], color=c2 )
    
   

    #leg = plt.legend(loc=(1.1, 0.4))
    #leg.get_frame().set_linewidth(0.0)
    
def plot_perturbed_channel_reference_cell(fig, pos, make_title=True, mesh=False):
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.spines['bottom'].set_visible(False)
    sub.spines['left'].set_visible(False)
    sub.set_xlim(0,9.)
    sub.set_ylim(0,8.8)
    sub.set_xticks([])
    sub.set_yticks([])
    sub.fill_between([0,9.],[0,0],[8.8, 8.8],color='steelblue', label='Liquid phase', lw=0.)
    
    sub.fill_between([0., 9.],[0.,0.], [1.2, 1.2], color='darkgrey', label='Solid phase', lw=0.)
    sub.fill_between([4., 5.],[1.2,1.2], [5.2, 5.2], color='darkgrey', lw=0.)
    sub.fill_between([0., 2.],[3.6,3.6], [8.8, 8.8], color='darkgrey', lw=0.)
    sub.fill_between([7., 9.],[3.6,3.6], [8.8, 8.8], color='darkgrey', lw=0.)
    sub.fill_between([2., 7.],[7.6,7.6], [8.8, 8.8], color='darkgrey', lw=0.)
    
    if make_title == True:
        sub.set_title('Perturbed channel 2D', fontsize=8)
    
    if mesh == True:
        lw_mesh = 0.2
        col_mesh = 'k'
        for i in range(18):
            for j in range(22):
                x1 = 0.5 * i
                x2 = 0.5 * (i + 1)
                y1 = 0.4 * j
                y2 = 0.4 * (j + 1)
                sub.plot([x1, x2],[y1, y1], lw=lw_mesh, c=col_mesh)
                sub.plot([x1, x1],[y1, y2], lw=lw_mesh, c=col_mesh)
                sub.plot([x1, x2],[y2, y1], lw=lw_mesh, c=col_mesh)
        sub.plot([0.,9.],[8.8,8.8], lw=lw_mesh, c=col_mesh)
        sub.plot([9.,9.],[0.,8.8], lw=lw_mesh, c=col_mesh)
    
    #leg = plt.legend(loc=(1.1, 0.4))
    #leg.get_frame().set_linewidth(0.0)
    
def plot_vertical_filaments_reference_cell(fig, pos):
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.spines['bottom'].set_visible(False)
    sub.spines['left'].set_visible(False)
    sub.set_xlim(0,36.)
    sub.set_ylim(0,36.)
    sub.set_xticks([])
    sub.set_yticks([])
    sub.fill_between([0,36.],[0,0],[36., 36.],color='steelblue', label='Liquid phase')
    
    sub.fill_between([5., 31.],[5.,5.], [13., 13.], color='darkgrey', label='Solid phase')
    sub.fill_between([0., 13.],[23.,23.], [31., 31.], color='darkgrey')
    sub.fill_between([23., 36.],[23.,23.], [31., 31.], color='darkgrey')
    
    sub.set_title('Vertical filaments 2D', fontsize=8)
    
    #leg = plt.legend(loc=(1.1, 0.4))
    #leg.get_frame().set_linewidth(0.0)
    
def plot_circle_in_a_box_reference_cell(fig, pos):
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.spines['bottom'].set_visible(False)
    sub.spines['left'].set_visible(False)
    sub.set_xlim(0,36.)
    sub.set_ylim(0,36.)
    sub.set_xticks([])
    sub.set_yticks([])
    sub.fill_between([0,36.],[0,0],[36., 36.],color='steelblue', label='Liquid phase')
    
    x = np.linspace(9., 27., 1000, endpoint=True)
    y1 = np.sqrt(np.square(9.) - np.square(x-18.)) + 18.
    y2 = -np.sqrt(np.square(9.) - np.square(x-18.)) + 18.
    
    sub.fill_between(x, y1, y2, color='darkgrey', label='Solid phase')
    
    sub.set_title('Circle in a box 2D', fontsize=8)
    
    #leg = plt.legend(loc=(1.1, 0.4))
    #leg.get_frame().set_linewidth(0.0)
    
def plot_single_box_in_a_box_reference_cell(fig, pos):
    sub = fig.add_axes(pos)
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    sub.spines['bottom'].set_visible(False)
    sub.spines['left'].set_visible(False)
    sub.set_xlim(0,36.)
    sub.set_ylim(0,36.)
    sub.set_xticks([])
    sub.set_yticks([])
    sub.fill_between([0,36.],[0,0],[36., 36.],color='steelblue', label='Liquid phase')
    
    x = np.linspace(10., 26., 10, endpoint=True)
    y1 = [10.] * 10
    y2 = [26.] * 10
    
    sub.fill_between(x, y1, y2, color='darkgrey', label='Solid phase')
    
    sub.set_title('Box in a box 2D', fontsize=8)
    
    #leg = plt.legend(loc=(1.1, 0.4))
    #leg.get_frame().set_linewidth(0.0)

    
def plot_pourosity_vs_correction_tensors_cap(fig, pos1, pos2):
    
    
    
    def e11(theta, eps_l, eps_s):
        return eps_l*theta + eps_s  * (1.-theta)
        
    def e22(theta, eps_l, eps_s):
        return eps_l / (theta + eps_l / eps_s * (1.-theta))
    
    pourosity = np.linspace(0., 1., 1000)
    
    
    sub1 = fig.add_axes(pos1)
    sub1.spines['top'].set_visible(False)
    sub1.spines['right'].set_visible(False)
    sub1.set_xlabel(r'Porosity $\theta$')
    sub1.set_ylabel(r'($E_{11}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$) [a. u.]')
    sub1.set_xlim([-0.05, 1.05])   
    sub1.set_ylim([-0.05, 1.05])   
    
    
    
    
    sub1.plot(pourosity, e11(pourosity, 2., 1.)-1., '-', label='$E_{11}$')
    sub1.plot(pourosity, (e11(pourosity, 100., 1.)-1.)/99., '-', label='$E_{11}$')
    sub1.plot(pourosity, (e11(pourosity, .01, 1.)-1.)/(0.01-1.), '-', label='$E_{11}$')
    sub1.plot([0.,1.],[0.,1.], '--', color='black') 

    #sub1.set_xscale('log')
    #sub1.legend()
    #sub1.set_title('Diffusion')
    
    sub2 = fig.add_axes(pos2)
    sub2.spines['top'].set_visible(False)
    sub2.spines['right'].set_visible(False)
    sub2.set_xlabel(r'Porosity $\theta$')
    sub2.set_ylabel(r'($E_{22}-\varepsilon_s$)/($\varepsilon_l - \varepsilon_s$) [a. u.]')
    sub2.set_xlim([-0.05, 1.05])   
    sub2.set_ylim([-0.05, 1.05])   
    
    sub2.plot(pourosity, (e22(pourosity, 100., 1.)-1.)/99., '-', label=r'$\varepsilon_l / \varepsilon_l = 100$')
    sub2.plot(pourosity, (e22(pourosity, 10., 1.)-1.)/9., '-', label=r'$\varepsilon / \varepsilon_l = 10$')
    sub2.plot(pourosity, e22(pourosity, 2., 1.)-1., '-', label=r'$\varepsilon / \varepsilon_l = 2$')
    sub2.plot(pourosity, (e22(pourosity, .5, 1.)-1.)/(0.5-1.), '-', label=r'$\varepsilon / \varepsilon_l = 0.5$')
    sub2.plot(pourosity, (e22(pourosity, .1, 1.)-1.)/(0.1-1.), '-', label=r'$\varepsilon / \varepsilon_l = 0.1$')
    sub2.plot(pourosity, (e22(pourosity, .01, 1.)-1.)/(0.01-1.), '-', label=r'$\varepsilon / \varepsilon_l = 0.01$')
    sub2.plot([0.,1.],[0.,1.], '--', color='black') 
    
    sub2.legend(loc=(1.0, 0.2))
    
    fig.suptitle('Straight channels 2D')
    
def plot_2d_coordinate_system(fig, pos):
    sub = fig.add_axes(pos)
    sub.set_xlabel(r'$x_1$')
    sub.set_ylabel(r'$x_2$')
    sub.set_xticks([])
    sub.set_yticks([])
    sub.spines['top'].set_visible(False)
    sub.spines['right'].set_visible(False)
    
    
####### FUNCTION FOR THEORETICAL PREDICIONS OF CHANNEL REFERENCE CELLS

def P(eps_l, eps_s, s1):
    p11 =s1*eps_l + (1-s1)*eps_s
    p12 = 0.
    p21 = 0.
    p22 = 1./((eps_s*s1+(1-s1)*eps_l)/(eps_l*eps_s))
    return np.array([[p11, p12],[p21, p22]])

def R(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def R_inv(theta):
    return np.array([[np.cos(theta), np.sin(theta)],[-np.sin(theta), np.cos(theta)]])

def straight_channel_theory_permittivity(eps_l, eps_s, component=''):
    p = np.linspace(0., 1., 1000, endpoint=True)
    y_data = np.zeros(1000)
    for i in range(1000):
        P_i = P(eps_l, eps_s, p[i])
        if component == '11':
            y_data[i] = (P_i[0,0]-eps_s)/(eps_l-eps_s)
        if component == '12':
            y_data[i] = (P_i[0,1]-eps_s)/(eps_l-eps_s)
        if component == '21':
            y_data[i] = (P_i[1,0]-eps_s)/(eps_l-eps_s)
        if component == '22':
            y_data[i] = (P_i[1,1]-eps_s)/(eps_l-eps_s)
    return [p, y_data]

def straight_channel_theory_diffusion(component=''):
    p = np.linspace(0., 1., 1000, endpoint=True)
    y_data = np.zeros(1000)
    for i in range(1000):
        D_i = np.array([[p[i],0.],[0., 0.]]) 
        if component == '11':
            y_data[i] = D_i[0,0]
        if component == '12':
            y_data[i] = D_i[0,1]
        if component == '21':
            y_data[i] = D_i[1,0]
        if component == '22':
            y_data[i] = D_i[1,1]
    return [p, y_data]

def diagonal_channel_theory_permittivity(eps_l, eps_s, component=''):
    p = np.linspace(0., 1., 1000, endpoint=True)
    y_data = np.zeros(1000)
    for i in range(1000):
        P_i = P(eps_l, eps_s, p[i])
        R_P = np.tensordot(R(np.pi/4.), P_i, axes=1)
        R_P_R_inv = np.tensordot(R_P, R_inv(np.pi/4.), axes=1)
        if component == '11':
            y_data[i] = (R_P_R_inv[0,0]-eps_s)/(eps_l-eps_s)
        if component == '12':
            y_data[i] = (R_P_R_inv[0,1]-eps_s)/(eps_l-eps_s)
        if component == '21':
            y_data[i] = (R_P_R_inv[1,0]-eps_s)/(eps_l-eps_s)
        if component == '22':
            y_data[i] = (R_P_R_inv[1,1]-eps_s)/(eps_l-eps_s)
    return [p, y_data]
eps_l =100.
eps_s = 1.

def diagonal_channel_theory_diffusion(component=''):
    p = np.linspace(0., 1., 1000, endpoint=True)
    y_data = np.zeros(1000)
    for i in range(1000):
        D_i = np.array([[p[i],0.],[0., 0.]])
        R_D = np.tensordot(R(np.pi/4.), D_i, axes=1)
        R_D_R_inv = np.tensordot(R_D, R_inv(np.pi/4.), axes=1)
        if component == '11':
            y_data[i] = R_D_R_inv[0,0]
        if component == '12':
            y_data[i] = R_D_R_inv[0,1]
        if component == '21':
            y_data[i] = R_D_R_inv[1,0]
        if component == '22':
            y_data[i] = R_D_R_inv[1,1]
    return [p, y_data]
        
