"""Test FUSOR tools."""

import pytest

from fusor.exceptions import IDTranslationException
from fusor.tools import translate_identifier


def test_translate_identifier(fusor_instance):
    """Test that translate_identifier method works correctly."""
    expected = "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"
    identifier = translate_identifier(fusor_instance.seqrepo, "NM_152263.3")
    assert identifier == expected

    identifier = translate_identifier(fusor_instance.seqrepo, "refseq:NM_152263.3")
    assert identifier == expected

    # test non-default target
    identifier = translate_identifier(
        fusor_instance.seqrepo, "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT", "refseq"
    )
    assert identifier == "refseq:NM_152263.3"

    # test no namespace
    with pytest.raises(IDTranslationException):
        identifier = translate_identifier(fusor_instance.seqrepo, "152263.3")

    # test unrecognized namespace
    with pytest.raises(IDTranslationException):
        identifier = translate_identifier(
            fusor_instance.seqrepo, "fake_namespace:NM_152263.3"
        )
