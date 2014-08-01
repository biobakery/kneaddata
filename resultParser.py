import re 
import collections

class ReadCounter():
    def __init__(self, fileType, pattern=r'.*', combineName=lambda m: m):
        '''
        pattern: regex pattern to get the name of the read. Default matches
        everything. 
        files: list of files to count
        fieldNames: list of field names to store
        fileType: what kind of file are we counting? sam, fastq, etc
        combineName: a function for combining the regex match
        '''
        self.strPattern = pattern
        self.strFileType = fileType

        # counter for the organisms
        self.counter = collections.Counter()

        self.combineName = combineName

    def skip(self, strLine):
        return False

    def processLine(self, strLine):
        self.counter[strLine] += 1

    def processPairedEnd(self):
        pass

    def count(self, strFile, bSingleEnd):
        #if not (strFile.endswith(self.strFileType)):
        #   #raise Exception(strFile + " is not a " + self.strFileType + " file!")

        try:
            f = open(strFile, "r")
            for line in f:
                if self.skip(line):
                    pass
                else:
                    self.processLine(line)
        except IOError:
            print("Could not find file " + strFile)
            self.iTotalAligned = None
            self.counter = None

        if not bSingleEnd:
            self.processPairedEnd()

    def get(self):
        return self.counter


class SamCounter(ReadCounter):
    def __init__(self, pattern=r'.*', combineName=lambda m: m):
        ReadCounter.__init__(self, fileType="sam", pattern=pattern,
                combineName=combineName)
        self.iTotalAligned = 0

    def skip(self, strLine):
        return (strLine[0] == '@')

    def processLine(self, strLine):
        lstrSplitLine = strLine.split("\t")
        # alignment!
        if lstrSplitLine[2].strip() != '*':
            self.iTotalAligned += 1
            match = re.search(self.strPattern, lstrSplitLine[0])
            if match:
                strOrganism = self.combineName(match)
                self.counter[strOrganism] += 1

    def processPairedEnd(self):
        for strOrganism in (self.counter).iterkeys():
            if self.counter[strOrganism] % 2 != 0:
                print("Invalid sam file for paired end reads!")
                print("You may have some non-unique reads here!")
                print(self.counter)
                print(self.iTotalAligned)
            self.counter[strOrganism] /= 2

        self.iTotalAligned /= 2

    def get(self):
        counter = ReadCounter.get(self)
        return (counter, self.iTotalAligned)


class FastqCounter(ReadCounter):
    def __init__(self, pattern=r'.*', combineName=lambda m: m):
        ReadCounter.__init__(self, fileType="fastq", pattern=pattern,
                combineName=combineName)
        self.iLineCounter = 0
        self.iLineMod = 4
        self.iTotalReads = 0

    def skip(self, strLine):
        # only consider every 4th line
        bSkip = True
        if (self.iLineCounter % self.iLineMod) == 0:
            bSkip = False
            self.iLineCounter = 0
        self.iLineCounter += 1
        return bSkip

    def processLine(self, strLine):
        match = re.search(self.strPattern, strLine)
        if match:
            strOrganism = self.combineName(match)
            self.counter[strOrganism] += 1
            self.iTotalReads += 1
        else:
            print strLine

    def get(self):
        counter = ReadCounter.get(self)
        return (counter, self.iTotalReads)


class BMTOutCounter(FastqCounter):
    def __init__(self, pattern=r'.*', combineName=lambda m: m):
        FastqCounter.__init__(self, pattern=pattern,
                combineName=combineName)

    def skip(self, strLine):
        return False
