function xlsdelemptysheets(filename)
% deletes all empty sheets in an xls-file
%
% Description
%
% This function loops through all sheets and deletes those sheets that are
% empty. Can be used to clean a newly created xls-file after all results
% have been saved in it.
%
% References: Torsten Jacobsen, "delete standard excel sheet" #, 22 Apr 2007 8:31 am </WebX?14@@.ef549db>
%
% Input:
%
% filename: name of xls file
%
% Output:
%
% none
%
% See also XLSWRITE
%

%
%=====================================================================
% Version : 1.0
% Author : hnagel
% Date : 27/04/2007
% Tested : 02/05/2007 (DR)
%=====================================================================
%
    % Check whether the file exists
    if ~exist(filename,'file')
       error([filename ' does not exist !']);
    else
        % Check whether it is an Excel file
        typ = xlsfinfo(filename);
        if ~strcmp(typ,'Microsoft Excel Spreadsheet')
           error([filename ' not an Excel sheet !']);
        end
    end

    % If filename does not contain a "\" the name of the current path is
    % added to filename. The reason for this is that the full path is
    % required for the command "excelObj.workbooks.Open(filename)" to work
    % properly.
    if isempty(strfind(filename,'\'))
        filename = [cd '\' filename];
    end

    excelObj = actxserver('Excel.Application');
    excelWorkbook = excelObj.workbooks.Open(filename);
    worksheets = excelObj.sheets;
    sheetIdx = 1;
    sheetIdx2 = 1;
    numSheets = worksheets.Count;

    % Loop over all sheets
    while sheetIdx2 <= numSheets
        % Saves the current number of sheets in the workbook
        temp = worksheets.count;
        % Check whether the current worksheet is the last one. As there
        % always need to be at least one worksheet in an xls-file the last
        % sheet must not be deleted.
        if or(sheetIdx>1,numSheets-sheetIdx2>0)
            worksheets.Item(sheetIdx).Delete;
        end
        % Check whether the number of sheets has changed. If this is not
        % the case the counter "sheetIdx" is increased by one.
        if temp == worksheets.count;
            sheetIdx = sheetIdx + 1;
        end
        sheetIdx2 = sheetIdx2 + 1; % prevent endless loop...
    end
    excelWorkbook.Save;
    excelWorkbook.Close(false);
    excelObj.Quit;
    delete(excelObj);

end