#! /usr/bin/env python

from __future__ import print_function
import rivet, re

authors_emails = {}
for aname in rivet.AnalysisLoader.analysisNames():
    ana = rivet.AnalysisLoader.getAnalysis(aname)
    for au_em in ana.authors():
        au, em = None, None
        if "<" not in au_em:
            if "@" in au_em:
                au = au_em
                em = au_em
            else:
                au = au_em
        else:
            m = re.search("(.*)<(.*)>.*", au_em)
            if m:
                au, em = m.group(1).strip(), m.group(2).strip()
        if au or em:
            if em or au not in authors_emails:
                authors_emails[au] = em

for au, em in sorted(authors_emails.items()):
    print(u"{} <{}>".format(au, em).encode("utf-8"))
