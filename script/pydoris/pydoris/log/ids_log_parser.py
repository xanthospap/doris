#! /usr/bin/python

import datetime
import re
import sys

def parse_first_log_line(line):
    """ Expect a line of type:
      '    YELLOWKNIFE DORIS site description form'
    """
    m = re.match(r"\s* [a-zA-Z_-]+ DORIS site description form\s*", line)
    if not m or len(m.groups()) != 1:
        errmsg = 'ERROR. Failed to match expected pattern in fist line \'{:}\''.format(
            line.strip)
        raise RuntimeError(errmsg)
    return m.group(1)


IdsLogFileBlocks = {
    0: '0. Form',
    1: '1. Site location information',
    2: '2. DORIS antenna and reference point information',
    3: '3. DORIS beacons information',
    4: '4. ITRF coordinates and velocities of the current DORIS ref. point (YEMB)',
    5: '5. IERS colocation information',
    6: '6. Tide Gauge colocation information',
    7: '7. Local site ties',
    8: '8. Meteorological Instrumentation',
    9: '9. DORIS network contacts'}


def line_is_block_header(line):
    for k, v in IdsLogFileBlocks.items():
        if re.match(r"{:}".format(v), line.strip()):
            return k
    return None


class IdsLogFile:

    def __init__(self, filename):
        self.filename_ = filename
        self.istream_ = None

    def open_stream(self, pos=0):
        if not self.istream_:
            self.istream_ = open(self.filename_, 'r')
        self.istream_.seek(0)

    def go_to_block(self, block_nr):
        self.open_stream()
        strpattern = IdsLogFileBlocks[block_nr]
        line = self.istream_.readline()
        line_matched = False
        while line:
            if re.match(r"{:}".format(strpattern), line):
                line_matched = True
                break
            line = self.istream_.readline()
        return line_matched

    def parse_this_block(self, block_nr):
        ret = []
        cdic = {}

        line = self.istream_.readline()
        while line:
            if len(line):
                next_line_read = False

                if re.match(r"{:}.[0-9]+\s*".format(block_nr), line):
                    # new dictionary/sub-header
                    if cdic:
                        ret.append(cdic)
                        cdic = {}

                elif line.startswith("  "):
                    l = [l.strip() for l in line.split(":")]
                    ## we would expect that spliting on ':' would give us two
                    ## columns but alas, some string can contain more than one
                    ## ':', e.g.:
                    ## Notes                    : 2 Ghz PC offset : 2 mm north, 1 mm west
                    ## Protect against that!
                    if len(l) != 2:
                        if len(l) > 2:
                            l = [l.strip() for l in line.split(":",1)]
                        else:
                            errmsg = "ERROR. Failed to parse line \'{:}\'".format(
                                line.strip())
                            raise RuntimeError(errmsg)
                    cdic[l[0]] = l[1]
                    next_line_read = True
                    line = self.istream_.readline()
                    while len(line) > 10 and line.startswith("  ") and ":" not in line:
                        cdic[l[0]] += ' ' + line.strip()
                        line = self.istream_.readline()

                elif line_is_block_header(line):
                    ret.append(cdic)
                    break

                if not next_line_read:
                    line = self.istream_.readline()
            else:
                line = self.istream_.readline()

        return ret

    def parse_block_nr(self, block_nr):
        if not self.go_to_block(block_nr):
            errmsg = "ERROR. Failed to find block nr {:} in file: ".format(
                block_nr, self.filename_)
            raise RuntimeError(errmsg)
        return self.parse_this_block(block_nr)


def dump_log_antenna_info(logfn,fout=sys.stdout):
    log = IdsLogFile(logfn)
    lofd = log.parse_block_nr(2)
    # for every dictionary (aka sub-block) in block-2
    for subblock in lofd:
        name = subblock["Four character ID"]
        height_str = subblock["Height above ground mark"]
        m = re.match(r"([0-9]+.[0-9]+)?\s*m", height_str)
        ## Guard against the case where no explicit "m" (aka meters) character
        ## is present. Produce a warning though
        if not m:
            m = re.match(r"([0-9]+.[0-9]+)?\s*", height_str)
            if not m:
                errmsg = "ERROR. Failed to parse height from string \'{}\' (file {:})".format(height_str,logfn)
                raise RuntimeError(errmsg)
            else:
                errmsg = "WARNING. No meters (\'m\' character) explicitely recorded in height string: \'{:}\' (file {:}). Assuming meters!".format(height_str,logfn)
                print(errmsg,file=sys.stderr)
        if len(m.groups()) != 1:
            errmmsg = "ERROR. Failed to parse site height from field: [{:}]".format(
                height_str)
            raise RuntimeError(errmsg)
        height = float(m.group(1)) if m.group(1) else 0e0
        installed_at = datetime.datetime.strptime(
            subblock["Date installed (dd/mm/yy)"], "%d/%m/%Y")
        removed_at = datetime.datetime.strptime(
            subblock["Date removed (dd/mm/yy)"], "%d/%m/%Y") if subblock["Date removed (dd/mm/yy)"] != "" else datetime.datetime.max
        print("{:4s} {:} {:} {:9.6f}".format(name, installed_at.strftime(
            "%Y-%m-%d %H:%M:%S"), removed_at.strftime("%Y-%m-%d %H:%M:%S"), height), file=fout)


if __name__ == "__main__":
    import sys

    dump_log_antenna_info(sys.argv[1])
    #log = IdsLogFile(sys.argv[1])
    #lofd = log.parse_block_nr(2)

    # for d in lofd:
    #  print('Sub-Block:')
    #  for k,v in d.items():
    #    print('{:} : {:}'.format(k,v))
