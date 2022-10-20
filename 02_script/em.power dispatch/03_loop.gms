parameter
price_roll(year_focus,month,day,hour,n,*)
;

loop(daily_window,
    x_down = 1+(ord(daily_window)-3)*24 ;
    x_up   =   ord(daily_window)*24 ;
    x_focus= 1+(ord(daily_window)-2)*24  ;
    x_focus_up= (ord(daily_window)-1)*24 ;

solve ProKoMo using LP minimizing COSTS ;

         price_roll(year_focus,month,day,hour,n,daily_window)      =  sum( t$(map_YMDHT(year_focus,month,day,hour,t)and ord(t)>=x_focus and ord(t)<= x_focus_up), energy_balance.M(n,t))*scaling_objective*(-1) ;

)
;
