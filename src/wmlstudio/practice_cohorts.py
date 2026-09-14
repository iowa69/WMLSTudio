"""Pinned public practice cohorts, downloaded on request and never committed.

Sequence data is never vendored in this repository or shipped in the portable
ZIP. A cohort is a table of NCBI RefSeq assembly accessions whose exact bytes
are pinned by size, MD5 and SHA-256; the files are fetched only when a user
explicitly asks for them. No download function is called by analysis or
application startup.

These are teaching cohorts, not validation sets. No expected ST, cluster,
threshold or "correct answer" is shipped with them, and the organism labels are
the ones NCBI records for those assemblies - they were not independently
verified here. Agreement between a WMLSTudio result and any published
investigation of the same strains is not a validation of this software.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from .sequence import check_cancelled, file_sha256

COHORT_FORMAT_VERSION = 1
NCBI_BASE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/'
MAX_GENOME_BYTES = 16 * 1024 * 1024
MAX_COHORT_BYTES = 128 * 1024 * 1024
MAX_CHECKSUM_BYTES = 1024 * 1024
LICENSE_NOTICE = (
    'Assemblies are public NCBI RefSeq records retrieved from '
    'https://ftp.ncbi.nlm.nih.gov/genomes/all/. NCBI places no restrictions on the use or '
    'distribution of the data it hosts, but it does not hold their copyright and cannot grant '
    'rights on behalf of the depositing submitters; individual submitters may assert terms. '
    'Cite the assembly accession and the originating submitters, not WMLSTudio, when reusing '
    'these sequences.')
GENOME_FIELDS = ('accession', 'assembly_dir', 'genus', 'species', 'strain',
                 'gz_bytes', 'gz_md5', 'gz_sha256', 'fasta_sha256')
SHARED_CAVEATS = (
    'These are public complete genomes from NCBI RefSeq. They are a teaching cohort, not a '
    'validation set.',
    'No expected ST, cluster, threshold or "correct answer" is shipped with this cohort. '
    'WMLSTudio computes relatedness from the scheme you choose; agreement with any published '
    'investigation is not a validation of this software.',
    'Strain and organism names are the labels NCBI records for these assemblies. They were not '
    'independently verified here.',
    'These are finished reference assemblies. Draft assemblies from your own sequencing will '
    'have more missing loci and different comparison behaviour.',
    'A genus or species folder created on import is a filing decision, not a laboratory '
    'identification. Confirm the organism before any clinical interpretation.',
)
COHORTS = {
    'kpneumoniae-10': {
        'title': 'Ten Klebsiella pneumoniae genomes (single species)',
        'purpose': 'Demonstrate single-species filing, MLST/cgMLST comparison and the minimum '
                   'spanning forest on one organism.',
        'expect': 'Every isolate should be proposed as one organism, so the import creates a '
                  'single genus/species folder. Nothing here proves an outbreak.',
        'citation_context': (
            'Rows 1-8 are the NIH Clinical Center Klebsiella pneumoniae complete-genome series '
            'described by Snitkin et al., Sci Transl Med 2012 (doi:10.1126/scitranslmed.3004129) '
            'and Conlan et al., Sci Transl Med 2014 (doi:10.1126/scitranslmed.3009845). Those '
            'papers are background for the cohort, not an expected result for this software.',
            'Rows 9-10 (NTUH-K2044, HS11286) are widely used unrelated reference strains included '
            'as visible non-cluster controls.'),
        'caveats': SHARED_CAVEATS,
        'genomes': (
            ('GCF_000281535.2', 'GCF_000281535.2_ASM28153v2', 'Klebsiella', 'pneumoniae', 'KPNIH1',
             1703154, 'd7d37504e86e35125f4cdfdb53765e55',
             '09a855d07afd582bbe73ee5b4cdb690bb55570ad2e830af77c240b63881dd707',
             'f590122e2cfa20d011575f04526330e000cba2740e78d5a7f1178bf64bbc9005'),
            ('GCF_000281435.2', 'GCF_000281435.2_ASM28143v2', 'Klebsiella', 'pneumoniae', 'KPNIH10',
             1703255, '42967c67a7c47e7c80579aa26e8be819',
             '7925758ba0b1be067f9203118532f4c1f4b57889805fe615e898c95f8e9633ac',
             'f7cf6d14b965797683a9b1ff7899c2d68f297cec90e5f88559883ea156830e56'),
            ('GCF_000714675.1', 'GCF_000714675.1_ASM71467v1', 'Klebsiella', 'pneumoniae', 'KPNIH24',
             1691578, 'd7727e770546b22c34371e8e5f8cce50',
             '257333af07ea24a878fcb371348ccfefb65ec14d06e5b3393e8df9f2601031ac',
             'c5213cf278756d6e6671942ca1f3806f953399e05e97a19db5da03f8392bed19'),
            ('GCF_000695935.1', 'GCF_000695935.1_ASM69593v1', 'Klebsiella', 'pneumoniae', 'KPNIH27',
             1810433, 'd715fe3519ed5e502f3bade189283a8a',
             '0faa17c39024a482cc0df917a855a43a4b7f3f90c1a02f236fb0029e3b3c85cc',
             '434148238674b69d570636d834d4dc8aba3e5864656bfa55380e03c6e1d93380'),
            ('GCF_000784945.1', 'GCF_000784945.1_ASM78494v1', 'Klebsiella', 'pneumoniae', 'KPNIH29',
             1629522, '13c367fa49546640f1fafc38b86592a1',
             '727a455e45bf7415bf2e616bcb637ea366f92ebfbb6e55e3abc66a78fd5a262b',
             '43a697d691a78afb76ff78aa2fb881a7236bcc0c6d167883edde9edb60855e33'),
            ('GCF_000784985.1', 'GCF_000784985.1_ASM78498v1', 'Klebsiella', 'pneumoniae', 'KPNIH30',
             1619174, '60cf353aab72a6e8e33338fece2ddbbc',
             '586c9934c6280ce1bb18fd1e063b014db1c5bcd9adba8934fbdb2b1e529eb02e',
             '0ec30a85b9d42482b7e65d1e52eb02df22017e576b98ebc5812f4742d6509745'),
            ('GCF_000785005.1', 'GCF_000785005.1_ASM78500v1', 'Klebsiella', 'pneumoniae', 'KPNIH31',
             1615609, '3ddd5fb27409c96c625dca7134b1290f',
             '2019edde516426611c56cc68de490057e4fa76bc222367ff416ca8ee3e6eea88',
             'b677d73f5c022ea1138ae6bfff2068652e3af66a2a7f5b3afbb9ef028ca15f8e'),
            ('GCF_000775375.1', 'GCF_000775375.1_ASM77537v1', 'Klebsiella', 'pneumoniae', 'KPNIH33',
             1678577, 'be43e9fff3063e176fa831c62c5fb37d',
             '47cb34bb41e4ba45f5c798b60f336b2b96c0710eb06ecc1cb751749c536edcd7',
             'a7722d088c6769f39369f79ac4181ba1f22f6285af44c1484931e24b57f104ca'),
            ('GCF_000009885.1', 'GCF_000009885.1_ASM988v1', 'Klebsiella', 'pneumoniae',
             'NTUH-K2044', 1615856, '79bef516fb4908294eaaf10d85c708ab',
             '9f24071c61d11e3df2fbe8fd2e56306afc90f1a13dc74a5ef6c4c80f014296a3',
             'd22693791a42bf78293a56366d525fc0428ee29fbd307ab3638452fa1c034292'),
            ('GCF_000240185.1', 'GCF_000240185.1_ASM24018v2', 'Klebsiella', 'pneumoniae', 'HS11286',
             1679403, 'd5ca66f37a3e464d1fa922748cfb87ff',
             'a6f9f6e36a6a891fd285b62da22543182dbd4ec728460958907edbc62831855c',
             '9528cd6940c58fae4b812a4af1a4c57f0db66f00b5ba3d60e0bb14609f423b8f'),
        )},
    'mixed-genus-20': {
        'title': 'Twenty genomes across fifteen genera',
        'purpose': 'Demonstrate automatic genus/species proposals and the folders an import '
                   'creates from them, including the calls the app should refuse to make.',
        'expect': 'Most isolates should be proposed confidently, but three are deliberately hard: '
                  'Klebsiella variicola is a Klebsiella pneumoniae species-complex member, '
                  'Escherichia coli cannot be separated from Shigella by whole-genome ANI, and '
                  'Enterobacter cloacae and Serratia marcescens have no species reference in the '
                  'bundled panel. Those should surface as a lower confidence tier or as needs-'
                  'review, never as a confident species call.',
        'citation_context': (
            'All twenty are public NCBI RefSeq complete genomes and are widely used as species '
            'reference strains. None of them is one of the assemblies in the bundled Klebsiella '
            'characterization panel, so no comparison in this cohort is a self-match against its '
            'own reference.',),
        'caveats': SHARED_CAVEATS + (
            'Klebsiella variicola belongs to the Klebsiella pneumoniae species complex. Separating '
            'it from K. pneumoniae is a genuinely hard call; treat a complex-level answer as the '
            'honest one unless you have independent evidence.',),
        'genomes': (
            ('GCF_000240185.1', 'GCF_000240185.1_ASM24018v2', 'Klebsiella', 'pneumoniae', 'HS11286',
             1679403, 'd5ca66f37a3e464d1fa922748cfb87ff',
             'a6f9f6e36a6a891fd285b62da22543182dbd4ec728460958907edbc62831855c',
             '9528cd6940c58fae4b812a4af1a4c57f0db66f00b5ba3d60e0bb14609f423b8f'),
            ('GCF_020525545.1', 'GCF_020525545.1_ASM2052554v1', 'Klebsiella', 'variicola', 'F2R9T',
             1628224, '1e97f644f527970d24476af2002622e8',
             '9284e4d25a225794a312891b11875f9a274e2f17e30a845fa6c5b122721d81b5',
             '194e5e0e71641aa7f61a59bc7b2eb67f35345833a814fae82e77312907c99269'),
            ('GCF_000008865.2', 'GCF_000008865.2_ASM886v2', 'Escherichia', 'coli', 'Sakai',
             1663694, 'e7bf4b6a0fa4cbd17a2cac69125311e0',
             '17bec26d56f9a188fa0c8aa2609be5a11ed7a0874450c9e02438f877907eb061',
             '71c2e5c364293c9ba36fc2c7acbcaa75cd6884295fe06260ba198826a8b1ddd3'),
            ('GCF_000013425.1', 'GCF_000013425.1_ASM1342v1', 'Staphylococcus', 'aureus', 'NCTC8325',
             820675, 'b086eb1020e7df022afa545dc6d93297',
             'f30fbb37ecacc1c6baca951ada10065a39a29e019c81f664a5d1fca3d8818adb',
             '05df047f78234de456a72df5f938b7fdef98fc3af3bf8634a11d70fbd3dcd21d'),
            ('GCF_000013465.1', 'GCF_000013465.1_ASM1346v1', 'Staphylococcus', 'aureus',
             'USA300-FPR3757', 847833, 'd2d2844a6812fd29d4e453069b93bac7',
             'cae163c4167f098bed2c0a6425351c0e6923904636e96cd4e44d792b71d5f9d0',
             'd5fec7a60a33287ea755c98754d9ca477341eede7aa9ab514cc0088a7458dd35'),
            ('GCF_000174395.2', 'GCF_000174395.2_ASM17439v2', 'Enterococcus', 'faecium', 'DO',
             898421, 'bf808e56028b29eafe5493981d987b02',
             'acbc39e584ef21030d28e93fa48a74d5a2971051ecef71e692507223a10eab78',
             'caa5e434e6fe6a28c8f5f0d59d620813d03963fa2065f8784b0d0cfe748a3aaf'),
            ('GCF_000007785.1', 'GCF_000007785.1_ASM778v1', 'Enterococcus', 'faecalis', 'V583',
             988302, '6176c3dee4e5c6b268902ed47b6ba5d5',
             'dfe7679cde7d62e00b35b14319164f17c2aff218baef9f1739ef2d672e0fbb9e',
             '57090755aa78ae8e3b65ccdddfffb9b9db41186ff2ebbb8dc5c8161b41a5c7ab'),
            ('GCF_000015425.1', 'GCF_000015425.1_ASM1542v1', 'Acinetobacter', 'baumannii',
             'ATCC17978', 1186598, '1e4ff3b78d6d9d11d9032697d2cd3517',
             'ac7ff42fa7bd113b70de8ec91920168af66e21fcb8d9297995b74a5f54e352b2',
             'c0b8caddc75fc55b5b208d2f4c328fec454962f7d2325c7df2c66fb04bd47d7c'),
            ('GCF_000006765.1', 'GCF_000006765.1_ASM676v1', 'Pseudomonas', 'aeruginosa', 'PAO1',
             1787220, 'c859ec5a25506367506bce972788aca8',
             '9433cb94b3f21fa40e94d0d9454996fcca39da1d49074b9b22b6a459e63ca510',
             '3e9335fb76772d324537b0c4c8cc7935e8b467b2acb0fae0c6e6dbeb2f227117'),
            ('GCF_000006945.2', 'GCF_000006945.2_ASM694v2', 'Salmonella', 'enterica', 'LT2',
             1471540, '3d4db127af17ce499dfec8cad09d4323',
             '2996f647c223d4c5b3e0f1a0940025614b39cad5c0c1a00b73f264e3b718d59a',
             '196a3eff4e6e6ed597173cc15d9dbde4cd55d76ad4d1b046da34cd75629c65cd'),
            ('GCF_000006885.1', 'GCF_000006885.1_ASM688v1', 'Streptococcus', 'pneumoniae', 'TIGR4',
             635909, '0c8e580789626aeb2e0d62d1bb6ddf28',
             'ba7ba0aa1ebe0752b83cf8d806ec497138dd2fc8020e2e0130dc4efea6e307d7',
             'd1de7c2eccd81e86d5ab24b4976b5ed867cd96d18ce17b6835f3b36a487036b7'),
            ('GCF_000196055.1', 'GCF_000196055.1_ASM19605v1', 'Streptococcus', 'agalactiae',
             'NEM316', 649835, '5c3233f4fef1bd00b7e5cda18dd47451',
             'a0f3ec657754efc0180580e1cacdb5e4b91e378ada56e6da874e2b28cd7b85d7',
             '9eedf8ae266186aab12c45d0b958dcbafda7837a186f9317f0ef55b21fce7c05'),
            ('GCF_000006785.2', 'GCF_000006785.2_ASM678v2', 'Streptococcus', 'pyogenes', 'SF370',
             548040, 'ca2a7ef6a1eb1f11be37e71e610e5f1a',
             '45baac35282c51956b5205e57869edc863b4e92bab8997e46ccc22932847ccd2',
             '99a7c6e07c6a0baa2d3deca583c07070b9f97d17ab46b692232f29ac6bb678e5'),
            ('GCF_000196035.1', 'GCF_000196035.1_ASM19603v1', 'Listeria', 'monocytogenes', 'EGD-e',
             869413, '79cf8e93793c9022e5648f8ff39cc75c',
             '753acf024fb30c338394467a0387e9449000ea2b2983123a9442e47e6c42c28c',
             'a0357eeedf295dd62f961b6f2c870f4c4c5ef55aac2fae03bcac38fc8c816292'),
            ('GCF_000008805.1', 'GCF_000008805.1_ASM880v1', 'Neisseria', 'meningitidis', 'MC58',
             648695, '287e55cbffed6f2874c691ad6717019a',
             '5e5a6ff7a4a84df95a3040a1865652f7eb47b02b8dfc7d9f55b76200d168b557',
             '88f5c3330d359c491607bcc9bd60afff7d41f4fad92ead675695f0556d183747'),
            ('GCF_000027305.1', 'GCF_000027305.1_ASM2730v1', 'Haemophilus', 'influenzae',
             'Rd-KW20', 539356, '22a41c5a66b3a51ecbba60061f4e365c',
             'c357cd91f6ca530879643ad19dbd90a976c25b137d8267bfdf63a3b0f7d6adac',
             '283170833d46e423548612f971189464fe6c96ab5a869b2f6bdf32b0eb69e1be'),
            ('GCF_000009205.2', 'GCF_000009205.2_ASM920v2', 'Clostridioides', 'difficile', '630',
             1214632, 'd8e93eecccfaff3f42e66d3ef490c0e6',
             '7671d3e15cde1da2b5ac87342e52bf308347d615ea7c38a68fd68d6dcb0418c2',
             'a0b0b1ccdc12fb4ad31c380517c0734764648b0c10d8828d784e562ae775c9ea'),
            ('GCF_000009085.1', 'GCF_000009085.1_ASM908v1', 'Campylobacter', 'jejuni', 'NCTC11168',
             465720, '8ca10c4de92c0c3fb2be58c9f81e32e9',
             '1a0c7e733310de29a6c715cf5e6825727f9ba3c9718b504cd5f5af03591599d9',
             'c64ee6a6db2c6407f8cca40a7048d1139770cd3f907a8c8dfb5c269c12b1567e'),
            ('GCF_000025565.1', 'GCF_000025565.1_ASM2556v1', 'Enterobacter', 'cloacae',
             'ATCC13047', 1662748, 'dde58dab6635faf25dd93864416918ee',
             '3e9169aadc55b680d36a3b5460cdfad7e7bfd043f84bf06cc1e12cf581180752',
             '20b07caba5dd2f35325c14913a46d2621e46a0e7a662fdab177c00bdb28fc749'),
            ('GCF_000513215.1', 'GCF_000513215.1_DB11', 'Serratia', 'marcescens', 'Db11',
             1497163, '4d25478af6345a51dd714260db0a72fc',
             '6c629eafa5806280dbbf814472bc4b2addac64bb2dd355a5eef9d7eeafe26d0c',
             '8219fa486c26aee5172e26610dfab158ad8360cb5bfc8b54c1b41eacd5412a13'),
        )},
}


def cohort_names():
    """List the pinned cohorts a user may ask for."""
    return tuple(sorted(COHORTS))


def cohort(name):
    if name not in COHORTS:
        raise ValueError(f'Unknown practice cohort: {name}')
    return COHORTS[name]


def cohort_genomes(name):
    """Return the pinned table as named records; the literals stay the only truth."""
    return [dict(zip(GENOME_FIELDS, row, strict=True)) for row in cohort(name)['genomes']]


def cohort_digest(name):
    """Fingerprint the pinned definition so an installed copy cannot silently drift."""
    content = {'format_version': COHORT_FORMAT_VERSION, 'cohort': name,
               'genomes': [list(row) for row in cohort(name)['genomes']]}
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def cohort_bytes(name):
    return sum(entry['gz_bytes'] for entry in cohort_genomes(name))


def describe_cohorts():
    """Plain summary rows for a listing CLI or a download dialog."""
    return [{'name': name, 'title': COHORTS[name]['title'], 'purpose': COHORTS[name]['purpose'],
             'expect': COHORTS[name]['expect'], 'genomes': len(COHORTS[name]['genomes']),
             'download_bytes': cohort_bytes(name),
             'genera': sorted({entry['genus'] for entry in cohort_genomes(name)}),
             'caveats': list(COHORTS[name]['caveats'])}
            for name in cohort_names()]


def _assembly_prefix(entry):
    accession = entry['accession']
    source, _, digits = accession.partition('_')
    digits = digits.split('.')[0]
    if not source or len(digits) != 9 or not digits.isdigit():
        raise ValueError(f'Unsupported assembly accession: {accession}')
    return f'{source}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{entry["assembly_dir"]}/'


def genome_relative_path(entry):
    """Path of the pinned genomic FASTA below the NCBI genomes/all/ root."""
    return _assembly_prefix(entry) + entry['assembly_dir'] + '_genomic.fna.gz'


def checksums_relative_path(entry):
    return _assembly_prefix(entry) + 'md5checksums.txt'


def genome_url(entry):
    return NCBI_BASE + genome_relative_path(entry)


def _safe_token(value):
    token = re.sub(r'[^A-Za-z0-9._-]+', '-', str(value)).strip('.-')
    return token or 'unknown'


def genome_filename(entry):
    """Name the download so the organism is readable but the folders are still unmade."""
    parts = (entry['genus'], entry['species'], entry['strain'], entry['accession'])
    return '_'.join(_safe_token(part) for part in parts) + '.fna.gz'


def default_destination(data_root, name):
    """Where a cohort installs under the user's Data root; never inside the app resources."""
    cohort(name)
    return Path(data_root) / 'practice-cohorts' / name


def cohort_manifest_path(path):
    return Path(path) / 'manifest.json'


def _open_source(relative, *, local_root=None, timeout=60):
    if local_root is not None:
        original = Path(local_root) / relative
        if original.is_symlink() or not original.is_file():
            raise ValueError(f'Practice cohort source must be a regular file: {original}')
        return original.open('rb')
    request = urllib.request.Request(NCBI_BASE + relative,
                                     headers={'User-Agent': 'WMLSTudio-practice-cohort-download'})
    return urllib.request.urlopen(request, timeout=timeout)


def _read_text(relative, *, local_root=None, cancelled=None, limit=MAX_CHECKSUM_BYTES):
    check_cancelled(cancelled)
    payload = b''
    with _open_source(relative, local_root=local_root) as handle:
        while chunk := handle.read(64 * 1024):
            check_cancelled(cancelled)
            payload += chunk
            if len(payload) > limit:
                raise ValueError(f'Practice cohort checksum file exceeds its bound: {relative}')
    return payload.decode('utf-8', 'replace')


def published_md5(entry, *, local_root=None, cancelled=None):
    """Read the MD5 NCBI publishes today for this assembly's genomic FASTA."""
    wanted = entry['assembly_dir'] + '_genomic.fna.gz'
    for line in _read_text(checksums_relative_path(entry), local_root=local_root,
                           cancelled=cancelled).splitlines():
        digest, _, path = line.strip().partition(' ')
        path = path.strip().lstrip('./')
        if path == wanted:
            return digest.strip().lower()
    return ''


def _fetch_genome(entry, target, *, local_root=None, cancelled=None):
    """Download one pinned genome and refuse to keep bytes that are not the pinned bytes."""
    relative = genome_relative_path(entry)
    digest = hashlib.sha256()
    md5 = hashlib.md5()
    total = 0
    with _open_source(relative, local_root=local_root, timeout=180) as handle, \
            target.open('xb') as output:
        while chunk := handle.read(1024 * 1024):
            check_cancelled(cancelled)
            total += len(chunk)
            if total > MAX_GENOME_BYTES:
                raise ValueError(f'Practice genome exceeds the 16 MiB download bound: {relative}')
            digest.update(chunk)
            md5.update(chunk)
            output.write(chunk)
    if total != entry['gz_bytes'] or digest.hexdigest() != entry['gz_sha256']:
        raise ValueError(f'Practice genome does not match its pinned SHA-256: {entry["accession"]}'
                         f' (expected {entry["gz_bytes"]} bytes / {entry["gz_sha256"]},'
                         f' received {total} bytes / {digest.hexdigest()})')
    if md5.hexdigest() != entry['gz_md5']:
        raise ValueError(f'Practice genome does not match its pinned MD5: {entry["accession"]}')
    upstream = published_md5(entry, local_root=local_root, cancelled=cancelled)
    if upstream != entry['gz_md5']:
        raise ValueError(f'NCBI publishes a different MD5 for {entry["accession"]}: the pinned '
                         f'record is stale and was not accepted (published {upstream or "nothing"},'
                         f' pinned {entry["gz_md5"]}).')
    return {'bytes': total, 'sha256': digest.hexdigest(), 'md5': md5.hexdigest()}


def _cached_genome(entry, cache, *, local_root=None, cancelled=None):
    """Reuse a content-addressed copy so an accession shared by two cohorts is fetched once."""
    cached = cache / f'{entry["gz_sha256"]}.fna.gz'
    if cached.is_file() and not cached.is_symlink() and cached.stat().st_size == entry['gz_bytes'] \
            and file_sha256(cached, cancelled) == entry['gz_sha256']:
        return cached, False
    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.practice-genome-', dir=cache) as temporary:
        staged = Path(temporary) / 'genome.fna.gz'
        _fetch_genome(entry, staged, local_root=local_root, cancelled=cancelled)
        check_cancelled(cancelled)
        os.replace(staged, cached)
    return cached, True


def download_cohort(name, destination, *, cache_root=None, source_root=None, cancelled=None,
                    progress=None):
    """Install one pinned practice cohort as flat files, or report the one already installed.

    Nothing partial is ever published: every genome is verified against its pinned size, MD5 and
    SHA-256 and against the MD5 NCBI publishes today, the whole cohort is staged in a temporary
    directory and only a complete, re-verified cohort is moved into place. The files are written
    flat and unsorted on purpose - creating the organism folders is the demonstration.
    """
    check_cancelled(cancelled)
    spec = cohort(name)
    entries = cohort_genomes(name)
    total = sum(entry['gz_bytes'] for entry in entries)
    if total > MAX_COHORT_BYTES:
        raise ValueError(f'Practice cohort exceeds the {MAX_COHORT_BYTES // (1024 * 1024)} MiB '
                         f'download bound: {name}')
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    cache = Path(cache_root).resolve() if cache_root is not None else destination.parent / '_cache'
    digest = cohort_digest(name)
    if destination.exists():
        existing = verify_cohort(destination, cancelled=cancelled)
        if existing['cohort_digest'] != digest:
            raise ValueError('A different practice cohort is already installed at that path; '
                             'choose a new destination.')
        return {'path': str(destination), 'cohort': name, 'cohort_digest': digest,
                'genomes': len(entries), 'downloaded': 0, 'reused': len(entries),
                'total_bytes': total, 'installed': False}
    manifest = {'format_version': COHORT_FORMAT_VERSION, 'cohort': name, 'cohort_digest': digest,
                'title': spec['title'], 'purpose': spec['purpose'], 'expect': spec['expect'],
                'source': 'NCBI RefSeq', 'source_base': NCBI_BASE,
                'source_mode': 'local_mirror' if source_root else 'pinned_https',
                'license_notice': LICENSE_NOTICE, 'caveats': list(spec['caveats']),
                'citation_context': list(spec['citation_context']),
                'created_utc': datetime.now(UTC).isoformat(), 'genomes': [], 'total_bytes': 0}
    downloaded = 0
    with tempfile.TemporaryDirectory(prefix='.practice-cohort-', dir=destination.parent) as work:
        stage = Path(work) / 'cohort'
        stage.mkdir()
        for index, entry in enumerate(entries, 1):
            check_cancelled(cancelled)
            cached, fetched = _cached_genome(entry, cache, local_root=source_root,
                                             cancelled=cancelled)
            downloaded += 1 if fetched else 0
            filename = genome_filename(entry)
            shutil.copyfile(cached, stage / filename)
            manifest['genomes'].append({
                'accession': entry['accession'], 'assembly_dir': entry['assembly_dir'],
                'genus': entry['genus'], 'species': entry['species'], 'strain': entry['strain'],
                'file': filename, 'bytes': entry['gz_bytes'], 'md5': entry['gz_md5'],
                'sha256': entry['gz_sha256'], 'fasta_sha256': entry['fasta_sha256'],
                'source_url': genome_url(entry),
                'checksum_basis': 'pinned_sha256_and_published_md5' if fetched
                                  else 'pinned_sha256_of_cached_copy'})
            if progress:
                progress(index, len(entries), f'Verified practice genome {entry["accession"]}')
        manifest['genomes'].sort(key=lambda item: item['file'])
        manifest['total_bytes'] = sum(item['bytes'] for item in manifest['genomes'])
        cohort_manifest_path(stage).write_text(json.dumps(manifest, indent=2) + '\n',
                                               encoding='utf-8')
        verify_cohort(stage, cancelled=cancelled)
        check_cancelled(cancelled)
        os.replace(stage, destination)
    return {'path': str(destination), 'cohort': name, 'cohort_digest': digest,
            'genomes': len(entries), 'downloaded': downloaded,
            'reused': len(entries) - downloaded, 'total_bytes': manifest['total_bytes'],
            'installed': True}


def verify_cohort(path, *, cancelled=None):
    """Re-hash an installed cohort against the pinned literals; raise on any disagreement."""
    path = Path(path).resolve()
    manifest_path = cohort_manifest_path(path)
    if manifest_path.stat().st_size > MAX_CHECKSUM_BYTES:
        raise ValueError('Practice cohort manifest exceeds the 1 MiB limit.')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    name = manifest.get('cohort')
    if manifest.get('format_version') != COHORT_FORMAT_VERSION or name not in COHORTS:
        raise ValueError('Practice cohort manifest is not a supported cohort record.')
    if manifest.get('cohort_digest') != cohort_digest(name):
        raise ValueError('Practice cohort manifest does not match the pinned cohort definition.')
    pinned = {entry['accession']: entry for entry in cohort_genomes(name)}
    seen = set()
    for item in manifest['genomes']:
        check_cancelled(cancelled)
        relative = item['file']
        target = path / relative
        if relative in seen or Path(relative).name != relative or target.is_symlink() \
                or not target.resolve().is_relative_to(path):
            raise ValueError('Unsafe or duplicated practice cohort path.')
        seen.add(relative)
        entry = pinned.get(item['accession'])
        if entry is None or item['sha256'] != entry['gz_sha256'] \
                or item['bytes'] != entry['gz_bytes'] or relative != genome_filename(entry):
            raise ValueError('Practice cohort manifest disagrees with the pinned accession: '
                             f'{item["accession"]}')
        if not target.is_file() or target.stat().st_size != entry['gz_bytes'] \
                or file_sha256(target, cancelled) != entry['gz_sha256']:
            raise ValueError(f'Practice cohort genome was changed or is missing: {relative}')
    if len(seen) != len(pinned):
        raise ValueError(f'Practice cohort is incomplete: {len(seen)} of {len(pinned)} genomes.')
    return manifest


def verify_pins(name, *, source_root=None, cancelled=None, progress=None):
    """Compare the pinned MD5s against NCBI's published checksums without downloading genomes."""
    entries = cohort_genomes(name)
    report = {'cohort': name, 'checked': 0, 'matched': [], 'mismatched': [], 'missing': []}
    for index, entry in enumerate(entries, 1):
        check_cancelled(cancelled)
        upstream = published_md5(entry, local_root=source_root, cancelled=cancelled)
        report['checked'] += 1
        if not upstream:
            report['missing'].append({'accession': entry['accession'],
                                      'url': NCBI_BASE + checksums_relative_path(entry)})
        elif upstream != entry['gz_md5']:
            report['mismatched'].append({'accession': entry['accession'],
                                         'pinned_md5': entry['gz_md5'], 'published_md5': upstream})
        else:
            report['matched'].append(entry['accession'])
        if progress:
            progress(index, len(entries), f'Checked published checksum {entry["accession"]}')
    report['stale'] = bool(report['mismatched'] or report['missing'])
    return report
