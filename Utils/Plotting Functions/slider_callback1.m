
function slider_callback1(src,eventdata,arg1)
val = get(src,'Value')*3;
set(arg1,'Position',[0 -val 1 4])
end

