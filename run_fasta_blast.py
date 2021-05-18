import os, sys, time, urllib, json
from xmltramp2 import xmltramp
import urllib.request as urllib2
import temp_search


class runFastaBlast:
    file = open('info')
    infoFile = json.load(file)
    path = infoFile['Path']
    url = 'https://www.ebi.ac.uk/Tools/services/rest/fasta'
    jobId = ''

    def __init__(self):
        pass

    @classmethod
    def runfasta(cls, params):
        requestUrl = cls.url + '/run/'
        # for homology
        requestData = urllib.parse.urlencode(params)
        try:
            req = urllib2.Request(requestUrl, None)
            reqH = urllib2.urlopen(req, requestData.encode(encoding='utf_8', errors='strict'))
            cls.jobId = str(reqH.read(), 'utf-8')
            reqH.close()
        except urllib2.HTTPError as ex:
            print(ex.read(), file=sys.stderr)
            raise
        return cls.jobId

    @classmethod
    def getResult(cls, name):
        # Check status and wait if necessary
        cls.clientPoll(cls.jobId)
        # Get available result types
        resultTypes = cls.serviceGetResultTypes(cls.jobId)
        for resultType in resultTypes:
            try:
                if str(resultType['fileSuffix']).strip() == 'xml':
                    # Derive the filename for the result
                    filename = name + '.blast'
                    # Get the result
                    result = cls.serviceGetResult(cls.jobId, str(resultType['identifier']))
                    if (str(resultType['mediaType']) == "image/png" or str(resultType['mediaType']) == "image/jpeg"):
                        fmode = 'wb'
                    else:
                        fmode = 'w'
                    # Write a result file
                    fh = open(os.path.join(cls.path, filename), fmode)
                    fh.write(result)
                    fh.close()
            except:
                pass

    @classmethod
    def restrequest(cls, url):
        try:
            req = urllib2.Request(url)
            # Make the request (HTTP GET).
            reqH = urllib2.urlopen(req)
            resp = reqH.read()
            contenttype = reqH.getheader("Content-Type")
            if (len(resp) > 0 and contenttype != "image/png;charset=UTF-8"
                and contenttype != "image/jpeg;charset=UTF-8"
                and contenttype != "application/gzip;charset=UTF-8"):
                result = str(resp, 'utf-8')
            else:
                result = resp
            reqH.close()
        # Errors are indicated by HTTP status codes.
        except urllib2.HTTPError as ex:
            # Trap exception and output the document to get error message.
            print(ex.read(), file=sys.stderr)
            raise
        return result

    @classmethod
    def clientPoll(cls, jobId):
        result = 'PENDING'
        while result == 'RUNNING' or result == 'PENDING':
            result = cls.serviceGetStatus(jobId)
            if result == 'RUNNING' or result == 'PENDING':
                time.sleep(10)

    @classmethod
    def serviceGetStatus(cls, jobId):
        requestUrl = cls.url + '/status/' + jobId
        status = cls.restrequest(requestUrl)
        return status

    @classmethod
    def serviceGetResultTypes(cls, jobId):
        requestUrl = cls.url + '/resulttypes/' + jobId
        xmlDoc = cls.restrequest(requestUrl)
        doc = xmltramp.parse(xmlDoc)
        return doc['type':]

    @classmethod
    def serviceGetResult(cls, jobId, type_):
        requestUrl = cls.url + '/result/' + jobId + '/' + type_
        result = cls.restrequest(requestUrl)
        return result

    # ___check previous job___

    @classmethod
    def serviceGetStatus(cls, jobId):
        requestUrl = cls.url + '/status/' + jobId
        status = cls.restrequest(requestUrl)
        return status

    # Print the status of a job
    @classmethod
    def printGetStatus(cls, jobId):
        status = cls.serviceGetStatus(jobId)
        print('Status:', status)
        if status == 'FINISHED':
            getInput = input('Get results (Y/N)? ')
            while True:
                if getInput.lower() == 'y':
                    cls.getResult(jobId, jobId)
                    break
                elif getInput.lower() == 'n':
                    break
                else:
                    continue


'''root = Tk()
fastaSet(root)
root.mainloop()

a = runFastaBlast
a.runfasta()
a.getResult(a.jobId)'''
