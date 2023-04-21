# Mathew Fischbach, fisch872@umn.edu

import logging
import sys


def get_logger(result_dir_path, loglevel):
    """Provide basic logging controls and returns a logging object.

    logging levels:
        DEBUG - Detailed information, useful for diagnosing problems
        INFO - Confirmation that things are working as expected
        WARNING - An indication that something unexpected happened, or a problem in the near future. still works as expected
        ERROR - Due to a more serious problem, the software has not been able to perform some function.
        CRITICAL - A serious error, indicating that the program itself may be unable to continue running.

    only call code if the logging is set at a certain level:
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug('Message with %s, %s', expensive_func1(),
                                                expensive_func2())
    """

    log_path = result_dir_path / 'logging'
    log_path.mkdir(exist_ok=True)

    logging._srcfile = None
    logging.logThreads = False
    logging.logProcesses = False
    logging.logMultiprocessing = False
    logging.basicConfig(filename=log_path / 'logger.log',
                        filemode='w',
                        format='%(levelname)s : %(asctime)s : %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=loglevel,
                        )
    logging.captureWarnings(True)
    logger = logging.getLogger(__name__)
    logger.info(sys.argv)
    return logger
