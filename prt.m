function prt(debg, level, txt, num)
% Print text and number to screen if debug is enabled.
if(debg >= level)
    if(numel(num) == 1)
        disp([txt num2str(num)]);
    else
        disp(txt)
        disp(num)
    end
end