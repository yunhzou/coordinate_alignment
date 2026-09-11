"""Keep presentation shared, and prevent display fitting from changing geometry."""
from pathlib import Path
import re

import numpy as np

from rxn_core.viewers import ASSETS, align_product, stylesheet


def test_only_two_shared_viewer_styles():
    assert sorted(p.name for p in ASSETS.glob('*viewer.css')) == [
        'catalog_viewer.css', 'reaction_viewer.css']
    assert '#ffd700' in stylesheet('reaction')
    assert 'background:#fafafa' in stylesheet('reaction')


def test_report_generators_do_not_define_independent_skins():
    root = Path(__file__).resolve().parents[1]
    for folder in ('bench', 'tools'):
        for path in (root / folder).glob('*.py'):
            if path.name == 'sync_viewer_styles.py':
                continue  # Recognizes existing style tags; it does not define CSS.
            assert not re.search(r'<style(?:>|\s)', path.read_text()), path
    for name in ('elementary_comparison_viewer.html', 'missing_patterns_viewer.html'):
        assert not (root / 'bench' / name).exists()


def test_display_fit_is_a_proper_rigid_transform():
    r = np.array([[0., 0., 0.], [2., 0., 0.], [.2, 1.3, 0.], [.1, .3, 1.7]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    p = r @ rotation + [9., -4., 2.]
    aligned = np.asarray(align_product(r, p))
    np.testing.assert_allclose(aligned, r, atol=1e-12)
    mirrored = p * [-1, 1, 1]
    fitted = np.asarray(align_product(r, mirrored))
    np.testing.assert_allclose(
        np.linalg.norm(fitted[:, None] - fitted[None, :], axis=-1),
        np.linalg.norm(mirrored[:, None] - mirrored[None, :], axis=-1), atol=1e-12)
    assert np.sign(np.linalg.det(fitted[1:] - fitted[0])) == np.sign(np.linalg.det(mirrored[1:] - mirrored[0]))
