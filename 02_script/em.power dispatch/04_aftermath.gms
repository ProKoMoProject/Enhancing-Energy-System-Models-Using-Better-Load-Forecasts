parameter

price(year_focus,month,day,hour,n)
;

price(year_focus,month,day,hour,'DE')               = sum(daily_window, price_roll(year_focus,month,day,hour,'DE',daily_window)     );



EXECUTE_UNLOAD '%output_dir%%result%.gdx'    price
;


$onecho >out.tmp

         par=price                         rng=price!A1:AJ9999    rdim=4 cdim=1

$offecho

execute "XLSTALK -c    %output_dir%%result%.xlsx" ;

EXECUTE 'gdxxrw %output_dir%%result%.gdx o=%output_dir%%result%.xlsx EpsOut=0 @out.tmp'
;
