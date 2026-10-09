# -*- coding: utf-8 -*-
# MGF -> mzML (1.1.0) 最小转换器：保留 TITLE/charge/RT/PEPMASS，供 Sage 等只读
# mzML 的引擎使用。谱 ID 用 "scan=<N>"（N 为顺序号，与 _renumber_mgf 的映射配套
# 或直接以 TITLE 顺序对齐）。
import sys, struct, base64

src, dst = sys.argv[1], sys.argv[2]

def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))

specs = []
cur = None
for line in open(src, encoding="utf-8", errors="replace"):
    ls = line.strip()
    if ls == "BEGIN IONS":
        cur = {"peaks": [], "charge": 2, "rt": None, "pepmass": None, "title": ""}
    elif ls == "END IONS":
        if cur is not None:
            specs.append(cur)
        cur = None
    elif cur is not None:
        if ls.startswith("TITLE="):
            cur["title"] = ls[6:]
        elif ls.startswith("CHARGE="):
            cur["charge"] = ls[7:].rstrip("+").split()[0]
        elif ls.startswith("RTINSECONDS="):
            cur["rt"] = ls[12:]
        elif ls.startswith("PEPMASS="):
            cur["pepmass"] = ls[8:].split()[0]
        elif ls and ls[0].isdigit():
            parts = ls.split()
            if len(parts) >= 2:
                cur["peaks"].append((parts[0], parts[1]))

with open(dst, "w", encoding="utf-8", newline="\n") as f:
    f.write('<?xml version="1.0" encoding="utf-8"?>\n'
            '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
            'xsi:schemaLocation="http://psi.hupo.org/ms/mzml '
            'http://psi.hupo.org/schema/mzml/1.1.0/mzML1.1.0.xsd">\n'
            '<mzML version="1.1.0">\n'
            '<cvList count="2">'
            '<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>'
            '<cv id="UO" fullName="Unit Ontology" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>'
            '</cvList>\n'
            '<fileDescription><fileContent/></fileDescription>\n'
            '<run id="run1">\n<spectrumList count="%d">\n' % len(specs))
    for n, sp in enumerate(specs, 1):
        mz_arr = base64.b64encode(struct.pack("<%dd" % len(sp["peaks"]),
                            *[float(p[0]) for p in sp["peaks"]])).decode()
        it_arr = base64.b64encode(struct.pack("<%dd" % len(sp["peaks"]),
                            *[float(p[1]) for p in sp["peaks"]])).decode()
        rt = sp["rt"] or "0"
        pep = sp["pepmass"] or "0"
        f.write('<spectrum index="%d" id="scan=%d" defaultArrayLength="%d">\n'
                % (n - 1, n, len(sp["peaks"])))
        f.write('<cvParam cvRef="MS" accession="MS:1000579" name="MSn spectrum" value=""/>\n')
        f.write('<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>\n')
        f.write('<scanList count="1"><scan>\n')
        f.write('<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="%s" unitName="second"/>\n' % rt)
        f.write('<userParam name="title" value="%s" type="xsd:string"/>\n' % esc(sp["title"]))
        f.write('</scan></scanList>\n')
        f.write('<precursorList count="1"><precursor>\n'
                '<selectedIonList count="1"><selectedIon>\n')
        f.write('<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="%s"/>\n' % pep)
        f.write('<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="%s"/>\n' % sp["charge"])
        f.write('</selectedIon></selectedIonList></precursor></precursorList>\n')
        f.write('<binaryDataArrayList count="2">\n')
        f.write('<binaryDataArray encodedLength="0">'
                '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit integer" value=""/>'
                '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit integer" value=""/>'
                '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>'
                '<binary>%s</binary></binaryDataArray>\n' % mz_arr)
        f.write('<binaryDataArray encodedLength="0">'
                '<cvParam cvRef="MS" accession="MS:1000521" name="32-bit integer" value=""/>'
                '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>'
                '<binary>%s</binary></binaryDataArray>\n' % it_arr)
        f.write('</binaryDataArrayList>\n</spectrum>\n')
    f.write('</spectrumList>\n</run>\n</mzML>\n</indexedmzML>\n')
print("wrote", len(specs), "spectra to", dst)
