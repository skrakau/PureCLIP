#!/usr/bin/env python3

import argparse
import re
import datetime
import sys

parser = argparse.ArgumentParser()
parser.add_argument('path', metavar = "<repository root>", help = "Path to PureCLIP repository root")
args = parser.parse_args()

def seqan_argparse_version():
    with open(args.path + '/src/pureclip.cpp') as pcpp:
        for line in pcpp:
            if re.match('.*setVersion\(parser,\s*"[.\d]+".*', line):
                return(line.split('"')[1].rstrip())
    raise Exception('Cannot detect pureclip version configuration of SeqAn argument parser')

def readme_version():
    with open(args.path + '/pkg/README') as readme:
        for line in readme:
            if re.match('^PureCLIP \d+[.]\d+[.]\d+$', line):
                return(line.split(' ')[1].rstrip())
    raise Exception('Cannot detect pureclip version configuration of README')

def readme_year():
    with open(args.path + '/pkg/README') as readme:
        for line in readme:
            if re.match('^\([cC]\) 20\d\d .*$', line):
                return(int(line.split(' ')[1].rstrip()))
    raise Exception('Cannot detect copyright year from README')

def cpack_version():
    major = None
    minor = None
    patch = None
    try:
        with open(args.path + '/src/CMakeLists.txt') as cmakelists:
            for line in cmakelists:
                if re.match('^(set|SET|Set)\(CPACK_PACKAGE_VERSION_MAJOR', line):
                    major = int(line.split('"')[1].rstrip())
                elif re.match('^(set|SET|Set)\(CPACK_PACKAGE_VERSION_MINOR', line):
                    minor = int(line.split('"')[1].rstrip())
                elif re.match('^(set|SET|Set)\(CPACK_PACKAGE_VERSION_PATCH', line):
                    patch = int(line.split('"')[1].rstrip())
        if all([ x != None for x in [major, minor, patch] ]):
            return(str(major) + '.' + str(minor) + '.' + str(patch))
    except ValueError:
        raise Exception('Cannot detect pureclip version configuration of CPack (ValueError)')
    raise Exception('Cannot detect pureclip version configuration of CPack')

####################################################################################################


try:
    success = True
    # Version Check
    versions = {
            'pureclip.cpp': seqan_argparse_version(),
            'README' : readme_version(),
            'CMakeLists.txt' : cpack_version()
            }
    if all([ v == versions['pureclip.cpp'] for v in versions.values()]):
        print("[OK] All files report version " + versions['pureclip.cpp'])
    else:
        print("[ERR] Ambiguous versions reported. " + format(versions), file = sys.stderr)
        success = False

    # Copyright Year Check
    cur_year = datetime.datetime.now().year
    if cur_year == readme_year():
        print("[OK] README copyright note contains current year " + str(cur_year))
    else:
        print("[ERR] README copyright note does not contain current year " + str(cur_year) + " but " + str(readme_year()), file = sys.stderr)
        success = False

    if success:
        sys.exit(0)
    else:
        sys.exit(1)

except Exception as err:
    print("[ERR] Package Sanity Check failed:\n" + format(err))
    sys.exit(1)




