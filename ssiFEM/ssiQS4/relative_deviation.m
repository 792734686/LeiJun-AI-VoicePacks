function errvalue = relative_deviation(value,value0)
    errvalue = 100*abs(value-value0)/max(abs(value0));
end