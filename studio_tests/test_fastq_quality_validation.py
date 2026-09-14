import pytest

from wmlstudio.sequence import _INVALID_QUALITY, SequenceReader


@pytest.mark.parametrize("value", range(128))
def test_compiled_quality_validation_matches_every_ascii_code(value):
    text = chr(value)
    original = any(ord(char) < 33 or ord(char) > 126 for char in text)
    assert bool(_INVALID_QUALITY.search(text)) == original
    assert bool(_INVALID_QUALITY.search("I" * 151 + text + "!" * 151)) == original


def test_compiled_quality_validation_rejects_non_ascii_without_range_guessing():
    for text in ("é", "😀", "\u0080", "\uffff"):
        assert _INVALID_QUALITY.search(text)


def test_all_permitted_quality_codes_remain_valid_in_real_wrapped_fastq(tmp_path):
    source = tmp_path / "all_qualities.fastq"
    quality = "".join(chr(value) for value in range(33, 127))
    source.write_text("@one\n" + "A" * len(quality) + "\n+\n" + quality[:40] + "\n" + quality[40:] + "\n", encoding="ascii")
    with SequenceReader(source) as reader:
        records = list(reader)
        assert reader.complete
    assert len(records) == 1 and records[0].quality == quality
