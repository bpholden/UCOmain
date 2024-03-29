from __future__ import print_function

import getpass
import os
import smtplib
import stat

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

import APF

# Possible logging levels:
# 'Info' 'Notice' 'Warn' 'Error'

APF.logident("master")


def logpush(filename, keep=4):
    """Takes the filename, and gives it the extension filename.1
       If there is alread a filename.1, it will become filename.2
       This continues untill filename.4, which will be deleted if it exists."""
    # First check that the initial file exists
    try:
        open(filename, 'r')
    except IOError:
        apflog("logpush: File %s could not be located." % filename, level="warn", echo=True)
        return
    # Try removing the oldest logged file. (currently number 4)
    try:
        os.remove(filename + ".4")
    except OSError:
        pass
    # Propagate the log extension down the files.
    for i in reversed(range(1, keep)):
        try:
            os.rename(filename + '.' + str(i), filename + '.' + str(i+1))
        except OSError:
            pass

    # change permissions
    try:
        os.chmod(filename, stat.S_IRUSR|stat.S_IRGRP|stat.S_IROTH|stat.S_IWUSR)
    except OSError:
        pass

    # Rename the main file filename.1
    try:
        os.rename(filename, filename + '.1')
    except OSError:
        apflog("Error renaming the primary file to be logged: %s."% (filename), \
             level="Error", echo=True)



def apflog(msg, level='Notice', echo=True):
    """Wraps the APF.log function. Messages are logged as 'master'."""

    APF.log(str(msg), level=level, echo=echo)

    if level in ['error']:
        subject = "[APF] An Error has occured"
        to_val = ['holden@ucolick.org', 'jrees@ucolick.org']
        sendmail(subject, msg, to=to_val)
    if level in ['Crit', 'Alert', 'Emerg']:
        subject = "[APF] A Serious Error has occured"
        to_val = ['holden@ucolick.org', '8314211210@txt.att.net', 'jrees@ucolick.org']
        sendmail(subject, msg, to=to_val)


def sendmail(subject, body, to=['holden@ucolick.org']):

    me = "APF <holden@ucolick.org>"

    msg = MIMEMultipart()
    msg["Subject"] = subject
    msg["From"] = me
    msg["To"] = ', '.join(to)
    msg.attach(MIMEText(body+"\n"))
    user_msg = "Script was being run by user %s" % getpass.getuser()
    msg.attach(MIMEText(user_msg))

    s = smtplib.SMTP('localhost')
    s.sendmail(me, to, msg.as_string())
    s.quit()

def main():
    body = "This is a test message. Error messages from the APF observe script will be sent with this function."
    subject = "[APF] Test"
    sendmail(subject, body, to=['holden@ucolick.org', '8314211210@txt.att.net'])

if __name__ == '__main__':
    # Send a test message
    main()
