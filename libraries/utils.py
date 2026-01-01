from config import Configuration
from datetime import datetime
import logging
import numpy as np
import os

class Utils:

    @staticmethod
    def create_directories(directory_list):
        '''This function checks for each directory in the list, and creates it, if it doesn't already exist.

        :parameter directory_list - The list of directories to create

        :return - Nothing is returned, but directories are created, if necessary
        '''

        for path in directory_list:
            # if the path does not exist then create it!
            if os.path.exists(path) is False:
                os.mkdir(path)
                Utils.log(path + ' created.', 'info')

    @staticmethod
    def get_all_files_per_date(path, file_ext):
        '''This function will return all the files for a directory split by date (typically for bias/dark/flat images).

        :parameter path - A string with the path of the file
        :parameter file_ext - The extension of the files, generally *.fits.

        :return files_no_ext - A list of files without their extension
        '''

        # get the files in the path with the given file extension
        files_no_ext = []
        for f in os.listdir(path):
            z = path + f + '/'
            try:
                for x in os.listdir(z):
                    if x.endswith(file_ext):
                        files_no_ext.append(z + x)
            except:
                Utils.log('No files found on ' + f, 'info')

        return files_no_ext

    @staticmethod
    def get_all_files_per_field(path, field, step, file_ext):
        '''This function will return all the files for a given TOROS field without their extension.

        :parameter path - A string with the path of the file
        :parameter field - The TOROS field
        :parameter step - Where in the process are you looking for files?
        :parameter file_ext - The type of file

        :return files_no_ext - A list of files without their extension
        :return uni_dte_dir - A list of the unique dates
        '''

        # get the files in the path with the given file extension
        files_no_ext = []
        sub_dir = os.listdir(path)
        dte_dir = []
        for f in sub_dir:
            if step == 'raw':
                z = path + f + '/' + field + '/' 
            else:
                z = path + f + '/' + field + '/'
            try:
                for x in os.listdir(z):
                    if ((x.split('_')[0] + '_' + x.split('_')[1] == field) |
                            (x.split('_')[0] == field)):
                        if x.endswith(file_ext):
                            files_no_ext.append(z + x)
            except:
                Utils.log('Field ' + field + ' not observed on ' + f, 'info')

        # get the unique dates with images
        for file in files_no_ext:
            dte_dir.append(file.split('/')[-3])
        uni_dte_dir = np.unique(dte_dir).tolist()

        return files_no_ext, uni_dte_dir

    @staticmethod
    def get_file_list(path, file_ext):
        '''This function will return the files in a given directory without their extension.

        :parameter path - A string with the path the file
        :parameter file_ext - The file type to make a list, generally *.fits

        :return file_list - A sorted list of files without the extension
        '''

        # get the files in the path with the given file extension
        file_list = [f for f in os.listdir(path) if f.endswith(file_ext)]

        # sort based on the number of the image
        file_list.sort(key=len)

        return file_list

    @staticmethod
    def is_time_between(begin_time, end_time):
        '''This function will check to see if the current time is between the provided times.

        :parameter begin_time - The star time for the check
        :parameter end_time - The end time for the check

        return - True or False
        '''

        # default to current UTC time
        check_time = datetime.now().time()

        if begin_time < end_time:
            return (check_time >= begin_time) & (check_time <= end_time)
        else:  # crosses midnight
            return (check_time >= begin_time) | (check_time <= end_time)

    @staticmethod
    def log(statement, level):
        '''This function logs all activity from the program to both screen and file.

        :parameter statement - A string which shows what needs to be logged
        :parameter level - The type of statement

        :return - Nothing is returned, but the log is updated and printed to the screen
        '''

        # create the logger
        logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', filename=Configuration.LOG_DIRECTORY + 'toros.log', filemode='a')
        logger = logging.getLogger()

        if not getattr(logger, 'handler_set', None):
            logger.setLevel(logging.DEBUG)

            # create console handler and set level to debug
            ch = logging.StreamHandler()
            ch.setLevel(logging.DEBUG)

            # create the formatter
            formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')

            # add formatter to ch
            ch.setFormatter(formatter)

            # add ch to logger
            logger.addHandler(ch)

            # 'set' Handler
            logger.handler_set = True

        if level == 'info':
            logger.info(statement)
        if level == 'debug':
            logger.debug(statement)
        if level == 'warning':
            logger.warning(statement)
        if level == 'error':
            logger.error(statement)
        if level == 'critical':
            logger.critical(statement)

    @staticmethod
    def write_txt(path, typ, line):
        '''This function will either write or append a text file, primarily used for results or logging.

        :parameter path: The location for the file
        :parameter typ: The write type either 'w', or 'a'
        :parameter line: The line you want to write to the file

        :return: Nothing is returned, but the file is written or appended.
        '''

        # open file
        file = open(path, typ)

        # write line to file
        file.write(line)

        # close file
        file.close()