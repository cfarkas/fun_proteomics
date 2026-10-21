import contextlib
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import fasta_download as db

FASTA = b'>sp|P1|PROTEIN_MOUSE OS=Mus musculus\nMPEPTIDER\n'


class DatabaseTests(unittest.TestCase):
    def test_local_target_and_gzip_auto_decoys(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
            source=Path('proteins.fasta.gz');source.write_bytes(gzip.compress(FASTA))
            cfg=db.normalize({'fasta':str(source), 'contaminants':'none'},tmp)
            db.build(cfg)
            rows=db.read_records('search.fasta')
            self.assertEqual(len(rows),2)
            self.assertTrue(rows[1][0].startswith('rev_'))
            self.assertEqual(cfg['sha256'],db.sha256(source))

    def test_prebuilt_is_kept_exactly_without_redoubling_decoys(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
            source=Path('prebuilt.fasta');source.write_bytes(FASTA+b'>rev_sp|P1|PROTEIN_MOUSE\nMEDITPEPR\n')
            cfg=db.normalize({'fasta':str(source),'contains_decoys':True,'contaminants':'none'},tmp)
            db.build(cfg)
            self.assertEqual(source.read_bytes(),Path('search.fasta').read_bytes())

    def test_partial_decoys_and_wrong_checksum_fail(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
            source=Path('bad.fasta');source.write_bytes(FASTA+b'>rev_sp|WRONG|WRONG\nMPEPTIDER\n')
            with self.assertRaisesRegex(ValueError,'SHA256'):
                db.normalize({'fasta':str(source),'sha256':'0'*64},tmp)
            cfg=db.normalize({'fasta':str(source),'contaminants':'none'},tmp)
            with self.assertRaisesRegex(ValueError,'exactly one'):
                db.build(cfg)

    def test_custom_contaminants_added_once_and_paired(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
            source=Path('proteins.fasta');source.write_bytes(FASTA)
            contaminants=Path('contaminants.fasta');contaminants.write_bytes(FASTA+b'>C2\nMPEPTIDEK\n')
            cfg=db.normalize({'fasta':str(source),'contaminants_fasta':str(contaminants)},tmp)
            db.build(cfg)
            rows=db.read_records('search.fasta')
            self.assertEqual(len(rows),4)
            self.assertEqual(len([h for h,s in rows if h.startswith('contam_')]),1)

    def test_source_selection_and_url_encoding(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'exactly one'):
                db.normalize({},tmp)
            with self.assertRaisesRegex(ValueError,'exactly one'):
                db.normalize({'taxon':9606,'proteome':'UP000005640'},tmp)
            cfg=db.normalize({'taxon':10090,'reviewed':False,'include_isoforms':True},tmp)
            self.assertIn('organism_id%3A10090',db.uniprot_url(cfg))
            self.assertNotIn('reviewed',db.uniprot_url(cfg))
            self.assertIn('includeIsoform=true',db.uniprot_url(cfg))
            with self.assertRaises(ValueError):db.require_http('file:///etc/passwd')

    def test_download_cache_offline_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings={'cache_dir':tmp,'retries':0}
            response=io.BytesIO(FASTA);response.headers={'X-UniProt-Release':'test'}
            with patch('urllib.request.urlopen',return_value=response) as request:
                path,info=db.download('https://example.org/ref.fasta',settings,validator=db.read_records)
                request.assert_called_once()
            with patch('urllib.request.urlopen',side_effect=AssertionError('Unexpected network')):
                cached,_=db.download('https://example.org/ref.fasta',{**settings,'offline':True},validator=db.read_records)
                self.assertEqual(path,cached)
            path.write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'SHA256'):
                db.download('https://example.org/ref.fasta',settings)
            with self.assertRaisesRegex(ValueError,'offline'):
                db.download('https://example.org/other.fasta',{**settings,'offline':True})

    def test_html_response_does_not_enter_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            response=io.BytesIO(b'<html>error</html>');response.headers={}
            with patch('urllib.request.urlopen',return_value=response), self.assertRaisesRegex(ValueError,'FASTA'):
                db.download('https://example.org/ref.fasta',{'cache_dir':tmp,'retries':0},validator=db.read_records)
            self.assertFalse(list(Path(tmp).glob('*.data')))
            self.assertFalse(list(Path(tmp).glob('*.partial')))

    def test_existing_contaminants_avoid_network(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
            source=Path('proteins.fasta');source.write_bytes(FASTA+b'>contam_C2\nMPEPTIDEK\n')
            cfg=db.normalize({'fasta':str(source)},tmp)
            with patch('urllib.request.urlopen',side_effect=AssertionError('Unexpected network')):
                db.build(cfg)
            self.assertEqual(len(db.read_records('search.fasta')),4)


if __name__=='__main__':unittest.main()
