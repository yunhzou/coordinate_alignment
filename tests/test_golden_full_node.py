from golden_full_node import aggregate


def test_cut_union_requires_all_negative_checks():
    no={'reference_recovery':'not_recovered'}
    assert aggregate([no],2)=='unknown'
    assert aggregate([no,no],2)=='not_recovered'
    assert aggregate([no,{'reference_recovery':'unknown'}],2)=='unknown'


def test_positive_cut_is_valid_in_partial_search():
    assert aggregate([{'reference_recovery':'recovered'}],50)=='recovered'
