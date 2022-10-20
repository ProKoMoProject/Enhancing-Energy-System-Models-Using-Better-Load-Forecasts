%% Data preprocessing of the TSO day-ahead load prediction for energy system models
% This script is used for data preprocessing of the TSO day-ahead load
% prediction, which is available on the ENTSO-E transparency platform.
% DATA:     TSO day-ahead load forecast and actual load (Excel file)
%           List of national holidays (Excel file)
% METHOD:   Using an decomposition model with seasonal components (weekday,
%           hour) and an SARIMA (1,0,1)(1,0,1)_24 
%           The whole model is estimated for every day-ahead forecast
%           iteratively and uses historical data from the last year for
%           estimation; the forecast is done as 24-hour-ahead forecast
% OUTPUT:   Improved load forecast

% created by: Mira Watermeyer; Department of Analytics and Statistics, Institute for Operations Research, Karlsruhe Institute for Technology
%             mira.watermeyer@kit.edu


clc
clear

%% Options: 
path = 'C:\Users\Mira\Documents\Data preprocessing\Paper\Code\';% path of load input
name_load = 'Load_DayAheadForecastAndActual_DL-LU_2016to2019';% name of file with load; file with 3 columns: time, forecast, actual
sheet_load = 'Sheet1'; % name of sheet with load data 
name_holiday = 'holiday_germany_2016to2019'; % list of holiday dates which have to be noted
sheet_holiday = 'holiday_nation'; % name of sheet with holiday dates 

start_time = datetime('01.01.16','Inputformat','dd.MM.yy','TimeZone','UTC'); %Timezone to capture time changes; Start date historical data

% Option to set for the user to decide how to treat holidays
holiday = 2; % 0: no holiday; 1: holiday as sunday (24&31.12. as saturday); 2: holiday as extra weekday (24&31.12. as saturday)

%%
try
%%                      *** Import holiday dates ***
opts = spreadsheetImportOptions("NumVariables", 1);
opts.Sheet = sheet_holiday;
opts.VariableNames = "Time";
opts.VariableTypes = "datetime";
opts = setvaropts(opts, "Time", "InputFormat", "");

% Import the data
holidaygermany = readtable([path, name_holiday, '.xlsx'], opts, "UseExcel", false);

clear opts


%% Loadprediction 

l_rw = calendarDuration(1,0,0);
Mdl_Time_Series = arima('ARLags',[1], 'MALags', [1], 'SARLags', [24], 'SMALags', [24]);
model = {'SARMA(1,1)(1,1)_24, rollingwindow 1y'}
 
    
%%                      *** Import load data ***
opts = spreadsheetImportOptions("NumVariables", 4);
opts.Sheet = sheet_load;

opts.VariableNames = ["Time", "Loadforecast", "Load"];
opts.VariableTypes = ["datetime", "double", "double"];
opts = setvaropts(opts, "Time", "InputFormat", "");

% Import the data
TimeLoadforecastLoadError = readtable([path, name_load, '.xlsx'], opts, "UseExcel", false);
TimeLoadforecastLoadError(1,:) = []; 
TimeLoadforecastLoadError.Error = TimeLoadforecastLoadError.Load - TimeLoadforecastLoadError.Loadforecast; 

clear opts

Data = timetable(TimeLoadforecastLoadError.Loadforecast,TimeLoadforecastLoadError.Load,TimeLoadforecastLoadError.Error,'TimeStep',hours(1),'StartTime',start_time, 'VariableNames',{'Loadforecast','Load','Error'});


%% MSE and MAE of forecasting error
mse_forecastingError = nanmean((Data.Error.^2));
mae_forecastingError = nanmean(abs(Data.Error));


%%                     *** Prepare data ***
% Find NaN entries
[row_nbr_missing_val_forecast, ] = find(isnan(Data(:,:).Loadforecast));
[row_nbr_missing_val_actual, ] = find(isnan(Data(:,:).Load));
date_missing_forecast = Data.Time(row_nbr_missing_val_forecast);
date_missing_actual = Data.Time(row_nbr_missing_val_actual);

% replace NaN: mean of load one week before and one week after  
for count = 1:length(row_nbr_missing_val_forecast)
    Data.Loadforecast(row_nbr_missing_val_forecast(count)) = nanmean([Data.Loadforecast(row_nbr_missing_val_forecast(count)-168); Data.Loadforecast(row_nbr_missing_val_forecast(count)+168)]);
end
for count = 1:length(row_nbr_missing_val_actual)
    Data.Load(row_nbr_missing_val_actual(count)) = nanmean([Data.Load(row_nbr_missing_val_actual(count)-168); Data.Load(row_nbr_missing_val_actual(count)+168)]);
end

Data.Error(isnan(Data(:,:).Error)) = Data.Load(isnan(Data(:,:).Error)) - Data.Loadforecast(isnan(Data(:,:).Error));




%%                    *** FORECAST ***
% DayAhead: 00:00 to 00:00 o'clock of following day; IntraDay: 12:00 to 00:00 o'clock of same day
numperiods_DayAhead = 24;
numperiods_IntraDay = 24;

end_date_training_data = start_time + l_rw - days(1) - hours(1); 
start_date_forecast_day_ahead = start_time + l_rw;
start_date_forecast_intraday = start_time + l_rw -days(1); 

[a ,forecast_start_date_day_ahead_rownbr] = max(Data.Time == start_date_forecast_day_ahead);
timevec_all_forecasts = Data(forecast_start_date_day_ahead_rownbr:end, :).Time;
forecast_day_ahead_summary = timetable(timevec_all_forecasts, zeros(size(timevec_all_forecasts)), zeros(size(timevec_all_forecasts)), zeros(size(timevec_all_forecasts)) , zeros(size(timevec_all_forecasts)), zeros(size(timevec_all_forecasts)), zeros(size(timevec_all_forecasts)), zeros(size(timevec_all_forecasts)), 'VariableNames',{'Forecast_Deseasonal','MSE_Forecast_Deseasonal', 'Seasonality', 'Errorforecast', 'IndexTrainDataStart', 'IndexTrainDataEnd','NewForecastError'} );

[a ,forecast_start_date_intra_day_rownbr] = max(Data.Time == start_date_forecast_intraday);
timevec_intra_day_forecast = Data(forecast_start_date_intra_day_rownbr:end, :).Time;
forecast_intra_day_summary = timetable(timevec_intra_day_forecast, zeros(size(timevec_intra_day_forecast)), zeros(size(timevec_intra_day_forecast)),zeros(size(timevec_intra_day_forecast)), 'VariableNames',{'Forecast_Deseasonal' ,'IndexTrainDataStart', 'IndexTrainDataEnd'} );

%capture parameter estimation in every loop iteration
param_estimates_summary_test = cell(2,5);
param_estimates_summary_test(1,:) = {'Time' , 'Intercept', 'AR', 'MA', 'Variance'};
param_estimates_summary = zeros(size(timevec_all_forecasts,1)/numperiods_DayAhead,4);

% capture seasonality in every loop iteration 
seasonality_loop = cell(1);
seasonality_loop(1,1) = {'Start TrainData'};
seasonality_loop(2,1) = {'End TrainData'};
seasonality_loop(3,1) = {'Start IntraDay'};
seasonality_loop(4,1) = {'Start DayAhead'};



%% BEGIN LOOP
j = days(0);
indicator = 0;
indicator_start_date = 0;
for i = 1:size(forecast_day_ahead_summary,1)/24 + 27
    % select training data for loop iteration
    end_date_training_data = end_date_training_data + j;
    
    % Case 1: end date exists in training data 
    if (ismember( end_date_training_data, Data.Time) == 1)
    % Case 3: end date in this iteration doesn't exist 
    elseif (ismember( end_date_training_data, Data.Time) == 0)
        memorized_end_date = end_date_training_data;
        end_date_training_data = end_date_training_data - caldays(1);
        indicator = 1;
        while (ismember(end_date_training_data, Data.Time) == 0)
            end_date_training_data = end_date_training_data - caldays(1);
            indicator = indicator + 1;
        end
    end 
    
    start_date_training_data = end_date_training_data - l_rw + hours(1); 

    % Logical check if start date exists
    % 1. Case: start date exists in data
    if (ismember( start_date_training_data, Data.Time) == 1)

    % 2. Case: start date doesn't exist for the first time
    elseif (ismember( start_date_training_data, Data.Time) == 0)
        % have to determine how many days to skip until start date is available
        k = 0;
        start_date_loop = start_date_training_data;
        while (ismember( start_date_loop, Data.Time) == 0)
            start_date_loop = start_date_loop + caldays(1);
            k = k + 1;
        end
        start_date_training_data = start_date_training_data + caldays(k);

    end

    [a ,forecast_start_date_rownbr] = max(Data.Time == start_date_training_data);
    [a ,forecast_end_date_rownbr] = max(Data.Time == end_date_training_data);
    training_Data_loop = Data(forecast_start_date_rownbr:forecast_end_date_rownbr,:);

    % forecast time vec
    start_date_forecast_day_ahead = start_date_forecast_day_ahead + j
    start_date_forecast_intraday = start_date_forecast_intraday + j
    forecast_vec_day_ahead = timetable(zeros(24,1),'TimeStep',hours(1),'StartTime',start_date_forecast_day_ahead, 'VariableNames',{'Seasonality'});
    forecast_vec_intraday = timetable(zeros(24,1),'TimeStep',hours(1),'StartTime',start_date_forecast_intraday, 'VariableNames',{'Seasonality'});

    is_member_logical = ismember(forecast_vec_day_ahead(:,:).Time,timevec_all_forecasts); % to check if forecast is even in ENTSOE data since we had missing values, don't predict them
    is_member_logical_intraday = ismember(forecast_vec_intraday(:,:).Time,timevec_intra_day_forecast);
    if is_member_logical == ones(24,1)

        % calculate seasonality of training data
        % seasonality depends on day of week and hour
        if holiday == 1 % incl. holidays
            [Data_seas_adj, forecast_vec_day_ahead] = myfun_adj_seas_incholiday(training_Data_loop, forecast_vec_day_ahead, holidaygermany);
        elseif holiday == 2 % holidays as extra weekday
            [Data_seas_adj, forecast_vec_day_ahead] = myfun_adj_seas_holiday_ex(training_Data_loop, forecast_vec_day_ahead, holidaygermany);
        else % only weekday and hour
            [Data_seas_adj, forecast_vec_day_ahead] = myfun_adj_seas(training_Data_loop, forecast_vec_day_ahead);
        end
        
        % save seasonality
        seasonality_loop(1,i+1) = {datestr(start_date_training_data)};
        seasonality_loop(2,i+1) = {datestr(end_date_training_data)};
        seasonality_loop(3,i+1) = {datestr(end_date_training_data + hours(13))};
        seasonality_loop(4,i+1) = {datestr(end_date_training_data + hours(25))};
        seasonality_loop(5,i+1) = {Data_seas_adj};
        
        % Use deseasonalized data and defined Model for forecast
        [EstMdl_Time_Series_loop, EstParamCov_loop, logL_loop, info_loop] = estimate(Mdl_Time_Series,Data_seas_adj(:,:).ForecastingError);
        [Y, YMSE] = forecast(EstMdl_Time_Series_loop, (numperiods_DayAhead + numperiods_IntraDay), Data_seas_adj(:,:).ForecastingError);
        
        if i == 1
            param_estimates_summary = (info_loop.X).';
        else
            param_estimates_summary(end+1, :) = (info_loop.X).';
        end
        is_member_logical_1 = ismember(timevec_all_forecasts, forecast_vec_day_ahead(:,:).Time);
        is_member_indices = find(is_member_logical_1);
        forecast_day_ahead_summary(is_member_indices,:).Seasonality = forecast_vec_day_ahead(:,:).Seasonality;
        forecast_day_ahead_summary(is_member_indices,:).Forecast_Deseasonal = Y(25:end);
        forecast_day_ahead_summary(is_member_indices,:).MSE_Forecast_Deseasonal = YMSE(25:end);
        forecast_day_ahead_summary(is_member_indices,:).IndexTrainDataStart = forecast_start_date_rownbr + zeros(24,1);
        forecast_day_ahead_summary(is_member_indices,:).IndexTrainDataEnd = forecast_end_date_rownbr + zeros(24,1);

        if is_member_logical_intraday == ones(24,1)
            is_member_logical_intraday_1 = ismember(timevec_intra_day_forecast, forecast_vec_intraday(:,:).Time);
            is_member_indices_intraday = find(is_member_logical_intraday_1);
            forecast_intra_day_summary(is_member_indices_intraday,:).Forecast_Deseasonal = Y(1:24);
            forecast_intra_day_summary(is_member_indices_intraday,:).IndexTrainDataStart = forecast_start_date_rownbr + zeros(24,1);
            forecast_intra_day_summary(is_member_indices_intraday,:).IndexTrainDataEnd = forecast_end_date_rownbr + zeros(24,1);
        end
    end
    j = caldays(1);

    if indicator ~= 0
        end_date_training_data = memorized_end_date;
        indicator = 0; % indicator is set to 0 again
    end
end

% create cell array for summary of parameter estimation
for i = 1:size(param_estimates_summary,2)
    param_estimates_summary_test(2,i+1) = {param_estimates_summary(:,i)};
end


% Getting all relevant data 
Data_forecastedTime = Data(find(start_time + l_rw == Data.Time):end,:);
forecast_day_ahead_summary.Errorforecast = forecast_day_ahead_summary.Forecast_Deseasonal + forecast_day_ahead_summary.Seasonality;
for t = 1:length(holidaygermany.Time)
    ind_feiertag = find((day(Data_forecastedTime.Time) == day(holidaygermany.Time(t))) & month(Data_forecastedTime.Time) == month(holidaygermany.Time(t)) & year(Data_forecastedTime.Time) == year(holidaygermany.Time(t)));
	forecast_day_ahead_summary.Errorforecast(ind_feiertag) = 0; 
end
forecast_day_ahead_summary.NewForecastError = Data(find(start_time + l_rw == Data.Time):end,:).Error - forecast_day_ahead_summary(:,:).Errorforecast;
forecast_day_ahead_summary.ImprLoadforecast = Data(find(start_time + l_rw == Data.Time):end,:).Load - forecast_day_ahead_summary(:,:).NewForecastError;

% Improvement in terms of MSE/MAE
mse_old = mean((Data(find(start_time + l_rw == Data.Time):end,:).Error.^2));
mse_new = mean(((Data(find(start_time + l_rw == Data.Time):end,:).Error - forecast_day_ahead_summary(:,:).Errorforecast).^2));

mae_old = mean(abs(Data(find(start_time + l_rw == Data.Time):end,:).Error));
mae_new = mean(abs(Data(find(start_time + l_rw == Data.Time):end,:).Error - forecast_day_ahead_summary(:,:).Errorforecast));

error_old_new_mat = [{''},{'old'},{'new'}; {'MSE'}, mse_old, mse_new; {'MAE'}, mae_old, mae_new]

% Save workspace of load perdiction 
save(['Loadprediction_Period', char(l_rw)]); 

catch
    warning(['Problem in estimation and forecast with used model, period size ', char(l_rw)]);
end








%% ------------------ FUNCTIONS ----------------------------

function [Data_seas_adj, forecast_vec] = myfun_adj_seas(Data, forecast_vec)

    code_vec = weekday(Data.Time) * 100 + hour(Data.Time);
    [code_unique,ia ,idx_code_unique] = unique(code_vec);
    mean_fe_per_code = accumarray(idx_code_unique, Data.Error,[],@mean);
    Data_seas_adj = zeros(size(Data.Time));
    for i = 1:size(code_unique,1)
         log = (code_vec == code_unique(i,:));
         val = Data.Error.*log;
         diff = (val - mean_fe_per_code(i,:)).*log;
         Data_seas_adj = Data_seas_adj + diff;
    end

    Data_seas_adj = timetable(Data.Time,Data_seas_adj, Data.Error - Data_seas_adj , 'VariableNames',{'ForecastingError', 'Seasonality'});

    code_forecast = weekday(forecast_vec.Time) * 100 + hour(forecast_vec.Time);
    seasonality_forecast = zeros(size(code_forecast));
    for i = 1:size(code_forecast,1)
        seasonality_forecast(i,:) = mean_fe_per_code(code_unique == code_forecast(i,:));
    end
    forecast_vec.Seasonality = seasonality_forecast;
end

function [Data_seas_adj, forecast_vec] = myfun_adj_seas_incholiday(Data, forecast_vec, holiday)
    
    weekday_data = weekday(Data.Time); 
    ind_2431 = find((day(Data.Time) == 24 & month(Data.Time) == 12)|(day(Data.Time) == 31 & month(Data.Time) == 12)); 
    ind_2431(weekday(Data.Time(ind_2431))==1) = []; 
    weekday_data(ind_2431) = 7; 
    
    for t = 1:length(holiday.Time)
        if ~(((day(holiday.Time(t))==24) && (month(holiday.Time(t))==12)) || ((day(holiday.Time(t))==31) && (month(holiday.Time(t))==12)))
            ind_feiertag = find((day(Data.Time) == day(holiday.Time(t))) & month(Data.Time) == month(holiday.Time(t)) & year(Data.Time) == year(holiday.Time(t)));
            weekday_data(ind_feiertag) = 1; 
        end
    end
    
    code_vec = weekday_data * 100 + hour(Data.Time);
    [code_unique,ia ,idx_code_unique] = unique(code_vec);
    mean_fe_per_code = accumarray(idx_code_unique, Data.Error,[],@mean);
    Data_seas_adj = zeros(size(Data.Time));
    for i = 1:size(code_unique,1)
         log = (code_vec == code_unique(i,:));
         val = Data.Error.*log;
         diff = (val - mean_fe_per_code(i,:)).*log;
         Data_seas_adj = Data_seas_adj + diff;
    end

    Data_seas_adj = timetable(Data.Time,Data_seas_adj, Data.Error - Data_seas_adj , 'VariableNames',{'ForecastingError', 'Seasonality'});

    weekday_data_forecast = weekday(forecast_vec.Time); 
    ind_2431 = find((day(forecast_vec.Time) == 24 & month(forecast_vec.Time) == 12)|(day(forecast_vec.Time) == 31 & month(forecast_vec.Time) == 12)); 
    ind_2431(weekday(forecast_vec.Time(ind_2431))==1) = []; 
    weekday_data_forecast(ind_2431) = 7; 
    
    for t = 1:length(holiday.Time)
        if ~(((day(holiday.Time(t))==24) && (month(holiday.Time(t))==12)) || ((day(holiday.Time(t))==31) && (month(holiday.Time(t))==12)))
            ind_feiertag = find((day(forecast_vec.Time) == day(holiday.Time(t))) & month(forecast_vec.Time) == month(holiday.Time(t)) & year(forecast_vec.Time) == year(holiday.Time(t)));
            weekday_data_forecast(ind_feiertag) = 1; 
        end
    end
    
    code_forecast = weekday_data_forecast * 100 + hour(forecast_vec.Time);
    seasonality_forecast = zeros(size(code_forecast));
    for i = 1:size(code_forecast,1)
        seasonality_forecast(i,:) = mean_fe_per_code(code_unique == code_forecast(i,:));
    end
    forecast_vec.Seasonality = seasonality_forecast;
end


function [Data_seas_adj, forecast_vec] = myfun_adj_seas_holiday_ex(Datainput, forecast_vec, holidaydate)
    
    weekday_data = weekday(Datainput.Time); 
    ind_2431 = find((day(Datainput.Time) == 24 & month(Datainput.Time) == 12)|(day(Datainput.Time) == 31 & month(Datainput.Time) == 12)); 
    ind_2431(weekday(Datainput.Time(ind_2431))==1) = 8; 
    weekday_data(ind_2431) = 8; 
    
    ind_feiertag_sum = 0; 
    for t = 1:length(holidaydate.Time)
        if ~(((day(holidaydate.Time(t))==24) && (month(holidaydate.Time(t))==12)) || ((day(holidaydate.Time(t))==31) && (month(holidaydate.Time(t))==12)))
            ind_feiertag = find((day(Datainput.Time) == day(holidaydate.Time(t))) & month(Datainput.Time) == month(holidaydate.Time(t)) & year(Datainput.Time) == year(holidaydate.Time(t)));
            weekday_data(ind_feiertag) = 8; 
            ind_feiertag_sum = ind_feiertag_sum + length(ind_feiertag); 
        end
    end
    
    code_vec = weekday_data * 100 + hour(Datainput.Time);
    [code_unique,ia ,idx_code_unique] = unique(code_vec);
    mean_fe_per_code = accumarray(idx_code_unique, Datainput.Error,[],@mean);
    Data_seas_adj = zeros(size(Datainput.Time));
    for i = 1:size(code_unique,1)
         log = (code_vec == code_unique(i,:));
         val = Datainput.Error.*log;
         diff = (val - mean_fe_per_code(i,:)).*log;
         Data_seas_adj = Data_seas_adj + diff;
    end

    Data_seas_adj = timetable(Datainput.Time,Data_seas_adj, Datainput.Error - Data_seas_adj , 'VariableNames',{'ForecastingError', 'Seasonality'});

    weekday_data_forecast = weekday(forecast_vec.Time); 
    ind_2431_forecast = find((day(forecast_vec.Time) == 24 & month(forecast_vec.Time) == 12)|(day(forecast_vec.Time) == 31 & month(forecast_vec.Time) == 12)); 
    ind_2431_forecast(weekday(forecast_vec.Time(ind_2431_forecast))==1) = 8; 
    weekday_data_forecast(ind_2431_forecast) = 8; 
    
    ind_feiertag_forecast_sum = 0; 
    for t = 1:length(holidaydate.Time)
        if ~(((day(holidaydate.Time(t))==24) && (month(holidaydate.Time(t))==12)) || ((day(holidaydate.Time(t))==31) && (month(holidaydate.Time(t))==12)))
            ind_feiertag_forecast = find((day(forecast_vec.Time) == day(holidaydate.Time(t))) & month(forecast_vec.Time) == month(holidaydate.Time(t)) & year(forecast_vec.Time) == year(holidaydate.Time(t)));
            weekday_data_forecast(ind_feiertag_forecast) = 8; 
            ind_feiertag_forecast_sum = ind_feiertag_forecast_sum + length(ind_feiertag_forecast); 
        end
    end
    
%     if (((isempty(ind_2431) == 1) && (ind_feiertag_sum == 0) && (isempty(find(weekday_data_forecast == 8))==0)) == 1)
%         if (isempty(ind_2431_forecast) == 0)
%             weekday_data_forecast(ind_2431_forecast) = 6; 
%         else
%             weekday_data_forecast(find(weekday_data_forecast == 8)) = 7; 
%         end
%             
%     end
    
    code_forecast = weekday_data_forecast * 100 + hour(forecast_vec.Time);
    seasonality_forecast = zeros(size(code_forecast));
    for i = 1:size(code_forecast,1)
        seasonality_forecast(i,:) = mean_fe_per_code(code_unique == code_forecast(i,:));
    end
    forecast_vec.Seasonality = seasonality_forecast;
end