
for ii in range(1,6):
    figure(ii)
    clf()

# Vector of half beam width [rad]
TH = array([2,10,20])*pi/180

# Angular values for evaluating beam intensity
th = linspace(-pi/2,pi/2,1000);

# MATLAB equiv.
#for ii = 1:length(TH):
#    theta_w = TH(ii);

for theta_w in TH:
    # Line array
    # Ratio of L/\lambda
    L_bar = 0.442/theta_w
    #th = theta_w
    x = pi*L_bar*sin(th)
    bl = (sin(x)/x)**2

    # Points for annotating polar plot
    th2 = linspace(-theta_w,theta_w,7)
    x = pi*L_bar*sin(th2)
    # Hack to avoid divide by zero
    x = x + 1e-3
    bl2 = (sin(x)/x)**2
    
    figure(1)
    plot(th/theta_w,10*log10(bl),
         label=r'%.1f deg, $L/\lambda$=%.1f'%(theta_w*180/pi,L_bar))

    figure(3)
    polar(th,10*log10(bl),
         label=r'%.1f deg, $L/\lambda$=%.1f'%(theta_w*180/pi,L_bar))

    figure(2)
    #plot(th/theta_w,10*log10(bl))
    plot(th/theta_w,bl,
         label=r'%.1f deg, $L/\lambda$=%.1f'%(theta_w*180/pi,L_bar))

    figure(4)
    #plot(th/theta_w,10*log10(bl))
    polar(th,bl,
         label=r'%.1f deg, $L/\lambda$=%.1f'%(theta_w*180/pi,L_bar))

    figure(5)
    polar(th/theta_w*30*pi/180.0,10*log10(bl),
         label=r'%.1f deg, $L/\lambda$=%.1f'%(theta_w*180/pi,L_bar))
    
figure(1)
axvline(-1,color='r',linestyle='--')
axvline(1,color='r',linestyle='--')
grid(True)
xlabel(r"Normalized angle, $\theta/\theta_{w}$ [n/a]")
ylabel(r"Source level [dB]")
ylim([-20,5])
legend()
tstr = r'$\theta_w$ = %.3f rad / %.1f deg, $L/\lambda$ = %.1f'%(theta_w,
                                                               theta_w*180/pi,
                                                               L_bar)
tstr = 'Line array beam pattern at different beam widths'
title(tstr)
xlim([-4,4])

figure(2)
axvline(-1,color='r',linestyle='--')
axvline(1,color='r',linestyle='--')
grid(True)
xlabel(r"Normalized angle, $\theta/\theta_{w}$ [n/a]")
ylabel(r"Source level [n/a]")
legend()
title(tstr)
xlim([-4,4])


figure(3)
xlim([-pi/2,pi/2])
xlim(45.0*pi/180*array([-1, 1]))
ylim([-20,3])
gca().set_theta_zero_location("N")
legend(loc='lower right',bbox_to_anchor=(0.2,0.1))
ylabel('dB')
yticks([-20, -15, -10, -5, -3, 0, 3])
xticks(arange(-40,50,10)*pi/180.0)

figure(4)
xlim([-pi/2,pi/2])
xlim(45.0*pi/180*array([-1, 1]))
#ylim([-20,5])
gca().set_theta_zero_location("N")
legend(loc='lower right',bbox_to_anchor=(0.2,0.1))
ylabel('Ratio [n/a]')
yticks([0, 0.25, 0.5, 0.75, 1.0])
xticks(arange(-40,50,10)*pi/180.0)


figure(5)
xlim([-pi/2,pi/2])
xlim(60.0*pi/180*array([-1, 1]))
ylim([-15,3])
gca().set_theta_zero_location("N")
#legend(loc='lower right',bbox_to_anchor=(0.2,0.1))
ylabel('dB')
yticks([-3,0])
xt = array([-30, 0, 30])*pi/180
xticks(xt,[r'$-\theta_w$',r'$0$',r'$\theta_w$'])

polar(th2/theta_w*30*pi/180.0,10*log10(bl2),'ro')
show()
