function o = format_tick(formatx, formaty)
x=str2num(get(gca,'XTickLabel'));
set(gca,'XTick',x)
set(gca,'XTickLabel',sprintf(formatx,x))
y=str2num(get(gca,'YTickLabel'));
set(gca,'YTick',y)
set(gca,'YTickLabel',sprintf(formaty,y))
end