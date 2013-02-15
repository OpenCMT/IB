#!/usr/bin/perl -w

# *************************************************************************
#
# This file is part of Qt Creator
#
# Copyright (c) 2012 Nokia Corporation and/or its subsidiary(-ies).
#
# Contact: Nokia Corporation (qt-info@nokia.com)
#
# GNU Lesser General Public License Usage
#
# This file may be used under the terms of the GNU Lesser General Public
# License version 2.1 as published by the Free Software Foundation and
# appearing in the file LICENSE.LGPL included in the packaging of this file.
# Please review the following information to ensure the GNU Lesser General
# Public License version 2.1 requirements will be met:
# http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
#
# In addition, as a special exception, Nokia gives you certain additional
# rights. These rights are described in the Nokia Qt LGPL Exception
# version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
#
# Other Usage
#
# Alternatively, this file may be used in accordance with the terms and
# conditions contained in a signed written agreement between you and Nokia.
#
# If you are unsure which license is appropriate for your use, please
# contact the sales department at http://qt.nokia.com/contact.
#
# *************************************************************************

use strict;
use Getopt::Long;
use IO::File;
use File::Copy;

my $USAGE=<<EOF;
Usage: generate.pl [--help] | [--recipe]


Custom wizard project generation example script.

EOF

my $argCount = scalar(@ARGV);
if ($argCount == 0
    || !GetOptions("help" => \$optHelp,
                   "recipe:s" => \$optRecipe)
    || $optHelp != 0) {
    print $USAGE;
    exit (1);
}

copy(\$optRecipe,"main.cpp") or die "Copy failed: $!";
