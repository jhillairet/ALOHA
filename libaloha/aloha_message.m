function form_str = aloha_message(str)
%  Format a string for information display
%  
%  form_str = aloha_message(str)
%  
%  form_str is of the form :
%  
%  '[ALOHA] (mmm.dd,yyyy HH:MM:SS) text of the string str'
%  
%  INPUTS
%  -------
%   - str [str] : input string message
%
%  OUTPUTS
%  -------
%   - form_str [str] : formatted output string message
%  
%  AUTHOR(S) : JH
%  
%  LAST UPDATE : 
%   - 17/06/2008 [creation]
%  
    msg_date = datestr(clock, 'mmm.dd,yyyy HH:MM:SS');
    
    form_str=['[ALOHA] (',msg_date,') ', str];