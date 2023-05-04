function settings = errormail(settings, errorstruct, i, parameterCell, evalfcn)

% depends on the sendmail command. Before this can be used, need to run the
% something like the following
%
% setpref('Internet','SMTP_Server','smtp.example.com');
% setpref('Internet','E_mail','your_email@example.com'');
% setpref('Internet','SMTP_Username','your_email@example.com');
% setpref('Internet','SMTP_Password','<replace me>'); 
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
%
% With whatever settings are apropriate for your email provider.
%

    persistent errorid
    
    if isempty(errorid) || ~strcmp(char(errorid), errorstruct.identifier)
        
        errorid = errorstruct.identifier;
        
        % Send email warning about this
        emailstr = sprintf(['This is an automated message from mcore\n\n', ...
            'The objective function has encountered the following error:\n\n%s\n\nWhile evaluating ind %d'], errorstruct.message, i);

        % convert the paramter evaluation furntion to a string so it can be
        % mentioned in the email report
        if isa(evalfcn, 'function_handle')
            evalfcn = func2str(evalfcn);
        else
            if ~ischar(evalfcn)
                warning('evalfcn was neither a string nor a function handle, mcore email not being sent.')
                return;
            end
        end

        if isempty (settings.ErrorFcnSettings.MailToAddress)
            warning('mcore objective fcn error email address in ErrorUserData.MailToAddress was empty. Email not sent, body was:\n\n%sn\n', ...
                emailstr);
        else
            try
                sendmail(settings.ErrorFcnSettings.MailToAddress, ...
                         sprintf('mcore: error in objective function %s', evalfcn), ...
                         emailstr);
            catch
                warning('mcore objective fcn error email send failure.');
            end
        end
    
    end

end