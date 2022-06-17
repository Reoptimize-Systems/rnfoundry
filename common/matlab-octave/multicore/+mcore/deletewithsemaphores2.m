function deletewithsemaphores2(fileList)
%DELETEWITHSEMAPHORES  Delete files using semaphors.
%   DELETEWITHSEMAPHORES(FILELIST) deletes the files in cell array FILELIST
%   using semaphores.
%
%		Example:
%		fileList = {'test1.txt', 'text2.txt'};
%		mcore.deletewithsemaphores(fileList);
%
%		Markus Buehren
%		Last modified 05.04.2009
%
%   See also SETFILESEMAPHORE, REMOVEFILESEMAPHORE.

    checkWaitTime = 0.1;
    nrOfAttempts  = 10;

    showWarnings  = 0;

    if ischar(fileList)
      fileList = {fileList};
    end

    for fileNr = 1:length(fileList)
        
      sem = mcore.setfilesemaphore(fileList{fileNr});
      
      for attemptNr = 1:nrOfAttempts
          
        throwError   = 0;
        
        showWarningsNow = showWarnings && (attemptNr == nrOfAttempts); % only in last attempt
        
        if mcore.mbdelete2(fileList{fileNr}, showWarningsNow, throwError)
          % deleting was successful
          break
        else
          % wait some time and check again
          pause(checkWaitTime);
        end
      end
      
      mcore.removefilesemaphore2(sem);
      
    end

end