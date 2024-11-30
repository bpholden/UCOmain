from __future__ import print_function

import getpass
import os
import smtplib
import stat
import datetime
import time

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


def timed_alert(subject, msg, to_pre, to_post):
    """Send an alert to the pre-defined list of people. If the current time is within the window
    of 10pm to 7am, the alert will be sent to the pre-defined list of people. Otherwise, it will be
    sent to the post-defined list of people.
    
    Useful to avoid waking people up in the middle of the night.
    """
    # Get the current time
    now = datetime.datetime.now()

    window_start = datetime.datetime(now.year, now.month, now.day, 23, 0, 0)
    window_end = datetime.datetime(now.year, now.month, now.day, 7, 0, 0)
    delta_day = datetime.timedelta(days=1)
    if now.hour < 7:
        window_start = window_start - delta_day
    else:
        window_end = window_end + delta_day

    # Check if the current time is within the window
    if now > window_start and now < window_end:
        to_val = to_post
    else:
        to_val = to_pre

    sendmail(subject, msg, to=to_val)

def apflog(msg, level='Notice', echo=True):
    """Wraps the APF.log function. Messages are logged as 'master'."""

    to_val = ['holden@ucolick.org', 'jrees@ucolick.org']
    alert_to_val = ['holden@ucolick.org', '8314211210@txt.att.net', 'jrees@ucolick.org']
    if level in ['error']:
        subject = "[APF] An Error has occured"
        sendmail(subject, msg, to=to_val)

    if level in ['Crit', 'Alert', 'Emerg']:
        subject = "[APF] A Serious Error has occured"
        sendmail(subject, msg, to=alert_to_val)

    if level in ['timed_alert']:
        subject = "[APF] A Serious Error has occured, but not so serious as to wake someone up."
        timed_alert(subject, msg, alert_to_val, to_val)
        level = 'Alert'

    APF.log(str(msg), level=level, echo=echo)

def sendmail(subject, body, to=['holden@ucolick.org']):
    """Send an email to the specified list of people."""

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
    body = "This is a test message."
    body += "Error messages from the APF observe script will be sent with this function."
    subject = "[APF] Test"
    to_val = ['holden@ucolick.org', '8314211210@txt.att.net']
    to_val_mod = to_val[0:1]
    sendmail(subject, body, to=to_val)
    time.sleep(3)
    timed_alert(subject, body, to_pre=to_val, to_post=to_val_mod)



if __name__ == '__main__':
    # Send a test message
    main()
