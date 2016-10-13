function [zout,pout] = pzplot_s(sys)
if iscell(sys)
    zout = sys{1};
    plot(zout,'rx','MarkerSize',12,'LineWidth',5)
    hold on
    pout = sys{2};
    plot(pout,'ro','MarkerSize',12,'MarkerFaceColor','r')
    sgrid%('ColorSpec',[0.1000 0.1000 0.1000])
    axis tight
    Xlimit = xlim;
    xlim([Xlimit(1) max(0,Xlimit(2))]);
    xlabel('Real axis($$sec^{-1}$$)')
    ylabel('Imaginary axis($$sec^{-1}$$)')
else
    try
      % Compute poles and zeros
      if issiso(sys)
         % Use data from ZPK representation in SISO case (for consistency
         % with ZPK(SYS) and PZPLOT/IOPZPLOT)
         [zout,pout] = zpkdata(sys,'v');
      else
         pout = pole(sys);
         zout = tzero(sys);
      end
    catch E
      throw(E)
    end
    % Call with graphical output
    try
        plot(pout,'ro','MarkerSize',12,'MarkerFaceColor','r')
        hold on
        plot(zout,'rx','MarkerSize',12,'LineWidth',5)
        sgrid%('ColorSpec',[0.1000 0.1000 0.1000])
        axis tight
        Xlimit = xlim;
        xlim([Xlimit(1) max(0,Xlimit(2))]);
        xlabel('Real axis($$sec^{-1}$$)')
        ylabel('Imaginary axis($$sec^{-1}$$)')
    catch E
      throw(E)
    end
end
end
