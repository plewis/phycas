# This utility provides a function (sendemail) that sends an email to a specified email address.
# One way this can be used is at the end of a Phycas Python script running on a remote server. After
# your analysis is finished, the following line will send a bodyless email to bunsen.honeydew@muppet.org
# with the sumject "Your Phycas run has finished":
# 
# from phycas.Utilities.enotify import sendemail
# sendemail(smtphost='smtp.muppet.org', toaddr='bunsen.honeydew@muppet.org', subject='Your Phycas run has finished')
#
# Notes: 
# - some STMP servers will reject an email if an invalid fromaddr (e.g. 'nobody') is specified
# - at this time, only one toaddr may be specified

import smtplib

def sendemail(smtphost, toaddr, subject, msgbody='', fromaddr='nobody'):
    msg = 'From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n%s\n' % (fromaddr, toaddr, subject, msgbody)
    server = smtplib.SMTP(smtphost)
    errmsg = None
    try:
        server.sendmail(fromaddr, toaddr, msg)
    except smtplib.SMTPRecipientsRefused:
        errmsg = 'SMTPRecipientsRefused exception'
    except smtplib.SMTPHeloError:
        errmsg = 'SMTPHeloError exception'
    except smtplib.SMTPSenderRefused:
        errmsg = 'SMTPSenderRefused exception'
    except smtplib.SMTPDataError:
        errmsg = 'SMTPDataError exception'

    if errmsg is None:
        import os
        if os.path.exists('smtperror.txt'):
            os.remove('smtperror.txt')
    else:
        open('smtperror.txt', 'w').write('From=%s\nTo=%s\nSubject=%s\nBody=%s\nSMTPhost=%s\nError=%s\n' % (fromaddr, toaddr, subject, msgbody, smtphost, errmsg))
        
    server.quit()
