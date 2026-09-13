"""Shared, presentation-only CSS injected across every Streamlit page.

Every rule here restyles how existing Streamlit elements render; it never
changes the text, order, or presence of any scientific content, and it is
safe to call from any page without affecting its data or computation.
"""

from __future__ import annotations

import streamlit as st

_GLOBAL_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

:root {
    --pee-brand: #2E7D46;
    --pee-brand-dark: #1B4D2C;
    --pee-brand-deep: #0F3320;
    --pee-brand-light: #5FBE7C;
    --pee-surface: #FFFFFF;
    --pee-surface-alt: #F1F6F2;
    --pee-border: #E1EAE3;
    --pee-text: #1F2A24;
    --pee-text-muted: #5B6B61;
}

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}

/* Remove the Streamlit Community Cloud "Deploy" affordance; this is a
finished, purpose-built app rather than a work-in-progress template. */
[data-testid="stAppDeployButton"] {
    display: none;
}

/* Sidebar navigation: subtle hover affordance only; Streamlit already
highlights the active page. */
[data-testid="stSidebarNavLink"] {
    border-radius: 8px;
    transition: background-color 0.15s ease;
}
[data-testid="stSidebarNavLink"]:hover {
    background-color: rgba(46, 125, 70, 0.10);
}

/* Alerts: replace Streamlit's saturated pastel fill with a quieter card
plus a coloured left accent, so scientific disclosures read as considered
notes rather than caution tape. */
[data-testid="stAlertContainer"] {
    border-radius: 10px;
    border: 1px solid var(--pee-border);
    border-left: 4px solid var(--pee-text-muted);
    background: var(--pee-surface-alt);
    padding: 0.9rem 1.1rem;
}
[data-testid="stAlertContainer"]:has([data-testid="stAlertContentInfo"]) {
    border-left-color: var(--pee-brand);
    background: #F1F8F3;
}
[data-testid="stAlertContainer"]:has([data-testid="stAlertContentSuccess"]) {
    border-left-color: var(--pee-brand);
    background: #F1F8F3;
}
[data-testid="stAlertContainer"]:has([data-testid="stAlertContentWarning"]) {
    border-left-color: #B7791F;
    background: #FBF6EC;
}
[data-testid="stAlertContainer"]:has([data-testid="stAlertContentError"]) {
    border-left-color: #B3261E;
    background: #FBEEED;
}
[data-testid="stAlertContainer"] p {
    color: var(--pee-text) !important;
}

/* Metrics as small cards with a brand-coloured top accent. */
[data-testid="stMetric"] {
    background: var(--pee-surface);
    border: 1px solid var(--pee-border);
    border-top: 3px solid var(--pee-brand);
    border-radius: 10px;
    padding: 0.85rem 1rem 0.65rem;
    box-shadow: 0 1px 2px rgba(15, 51, 32, 0.05);
}

/* Buttons: rounded, brand-weighted, with a light lift on hover. */
.stButton > button,
.stDownloadButton > button {
    border-radius: 999px;
    font-weight: 600;
    transition: transform 0.15s ease, box-shadow 0.15s ease;
}
.stButton > button:hover,
.stDownloadButton > button:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 10px rgba(15, 51, 32, 0.12);
}

/* Page-link grids (used on the home page): card-style rows instead of
plain text links. */
[data-testid="stPageLink"] a[data-testid="stPageLink-NavLink"] {
    background: var(--pee-surface);
    border: 1px solid var(--pee-border);
    border-radius: 10px;
    padding: 0.7rem 0.9rem;
    transition: transform 0.15s ease, box-shadow 0.15s ease, border-color 0.15s ease;
    box-shadow: 0 1px 2px rgba(15, 51, 32, 0.04);
    white-space: normal !important;
    height: auto !important;
}
[data-testid="stPageLink"] a[data-testid="stPageLink-NavLink"] span {
    white-space: normal !important;
    min-width: 0;
    word-break: break-word;
}
[data-testid="stPageLink"] a[data-testid="stPageLink-NavLink"]:hover {
    border-color: var(--pee-brand);
    box-shadow: 0 4px 10px rgba(15, 51, 32, 0.10);
    transform: translateY(-1px);
}

/* Bordered containers (st.container(border=True)): soften the default
grey border to match the rest of the palette. */
[data-testid="stVerticalBlockBorderWrapper"] {
    border-color: var(--pee-border) !important;
    border-radius: 12px !important;
}

.pee-eyebrow {
    display: inline-block;
    font-size: 0.72rem;
    font-weight: 700;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: var(--pee-brand);
    background: #E7F3EB;
    border-radius: 999px;
    padding: 0.3rem 0.8rem;
    margin-bottom: 0.6rem;
}

.pee-dark-panel {
    background: linear-gradient(135deg, var(--pee-brand-deep) 0%, var(--pee-brand-dark) 100%);
    border-radius: 20px;
    padding: 2rem 2.2rem;
    color: #EAF3EC;
    margin: 0.5rem 0 1.5rem;
}
.pee-dark-panel h2 {
    color: #FFFFFF;
    margin-top: 0;
}
.pee-dark-panel ol {
    list-style: none;
    counter-reset: pee-step;
    margin: 1.1rem 0 0;
    padding: 0;
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
    gap: 0.75rem;
}
.pee-dark-panel ol li {
    counter-increment: pee-step;
    position: relative;
    background: rgba(255, 255, 255, 0.07);
    border: 1px solid rgba(255, 255, 255, 0.14);
    border-radius: 12px;
    padding: 0.9rem 1rem 0.9rem 2.9rem;
    line-height: 1.45;
    font-size: 0.92rem;
}
.pee-dark-panel ol li::before {
    content: counter(pee-step);
    position: absolute;
    left: 0.8rem;
    top: 0.8rem;
    width: 1.5rem;
    height: 1.5rem;
    border-radius: 50%;
    background: var(--pee-brand-light);
    color: var(--pee-brand-deep);
    font-weight: 800;
    font-size: 0.8rem;
    display: flex;
    align-items: center;
    justify-content: center;
}
.pee-dark-panel strong {
    color: #FFFFFF;
}

.pee-chip-row {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
    margin: 0.75rem 0 1.25rem;
}
.pee-chip {
    background: var(--pee-surface-alt);
    border: 1px solid var(--pee-border);
    border-radius: 999px;
    padding: 0.4rem 0.9rem;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--pee-brand-dark);
}
</style>
"""


def inject_global_styles() -> None:
    """Inject the shared stylesheet. Call once near the top of every page."""

    st.html(_GLOBAL_CSS)


_HERO_ILLUSTRATION_SVG = (
    '<svg viewBox="0 0 480 360" xmlns="http://www.w3.org/2000/svg" role="img" '
    'aria-label="Stylised, illustrative preview of descriptive charts; not real data">'
    '<rect x="1" y="1" width="478" height="358" rx="20" fill="#FFFFFF" stroke="#E1EAE3"/>'
    '<path d="M1 21 V21 a20 20 0 0 1 20-20 h438 a20 20 0 0 1 20 20 v15 H1 Z" fill="#F1F6F2"/>'
    '<circle cx="20" cy="18" r="5" fill="#2E7D46"/>'
    '<circle cx="38" cy="18" r="5" fill="#5FBE7C"/>'
    '<circle cx="56" cy="18" r="5" fill="#B7D9C2"/>'
    '<line x1="24" y1="188" x2="230" y2="188" stroke="#DCE7E0" stroke-width="1.5"/>'
    '<line x1="24" y1="56" x2="24" y2="188" stroke="#DCE7E0" stroke-width="1.5"/>'
    '<circle cx="55" cy="140" r="7" fill="#2E7D46" opacity="0.85"/>'
    '<circle cx="72" cy="118" r="7" fill="#2E7D46" opacity="0.85"/>'
    '<circle cx="92" cy="150" r="7" fill="#2E7D46" opacity="0.85"/>'
    '<circle cx="62" cy="98" r="7" fill="#2E7D46" opacity="0.85"/>'
    '<circle cx="102" cy="130" r="7" fill="#2E7D46" opacity="0.85"/>'
    '<circle cx="172" cy="88" r="7" fill="#8FCB9F" opacity="0.9"/>'
    '<circle cx="196" cy="108" r="7" fill="#8FCB9F" opacity="0.9"/>'
    '<circle cx="160" cy="74" r="7" fill="#8FCB9F" opacity="0.9"/>'
    '<circle cx="207" cy="94" r="7" fill="#8FCB9F" opacity="0.9"/>'
    '<circle cx="186" cy="128" r="7" fill="#8FCB9F" opacity="0.9"/>'
    '<line x1="250" y1="188" x2="456" y2="188" stroke="#DCE7E0" stroke-width="1.5"/>'
    '<rect x="258" y="128" width="26" height="60" rx="3" fill="#2E7D46"/>'
    '<rect x="298" y="93" width="26" height="95" rx="3" fill="#4C9A63"/>'
    '<rect x="338" y="148" width="26" height="40" rx="3" fill="#2E7D46"/>'
    '<rect x="378" y="78" width="26" height="110" rx="3" fill="#5FBE7C"/>'
    '<rect x="418" y="118" width="26" height="70" rx="3" fill="#4C9A63"/>'
    '<g>'
    '<rect x="24" y="212" width="52" height="32" fill="#2E7D46" opacity="0.95"/>'
    '<rect x="78" y="212" width="52" height="32" fill="#4C9A63" opacity="0.8"/>'
    '<rect x="132" y="212" width="52" height="32" fill="#8FCB9F" opacity="0.65"/>'
    '<rect x="186" y="212" width="52" height="32" fill="#C9E3D1" opacity="0.55"/>'
    '<rect x="240" y="212" width="52" height="32" fill="#C9E3D1" opacity="0.45"/>'
    '<rect x="294" y="212" width="52" height="32" fill="#8FCB9F" opacity="0.55"/>'
    '<rect x="348" y="212" width="52" height="32" fill="#4C9A63" opacity="0.7"/>'
    '<rect x="402" y="212" width="54" height="32" fill="#2E7D46" opacity="0.9"/>'
    '<rect x="24" y="248" width="52" height="32" fill="#4C9A63" opacity="0.75"/>'
    '<rect x="78" y="248" width="52" height="32" fill="#2E7D46" opacity="0.95"/>'
    '<rect x="132" y="248" width="52" height="32" fill="#4C9A63" opacity="0.8"/>'
    '<rect x="186" y="248" width="52" height="32" fill="#8FCB9F" opacity="0.6"/>'
    '<rect x="240" y="248" width="52" height="32" fill="#C9E3D1" opacity="0.5"/>'
    '<rect x="294" y="248" width="52" height="32" fill="#C9E3D1" opacity="0.45"/>'
    '<rect x="348" y="248" width="52" height="32" fill="#8FCB9F" opacity="0.6"/>'
    '<rect x="402" y="248" width="54" height="32" fill="#4C9A63" opacity="0.75"/>'
    '<rect x="24" y="284" width="52" height="32" fill="#8FCB9F" opacity="0.6"/>'
    '<rect x="78" y="284" width="52" height="32" fill="#4C9A63" opacity="0.8"/>'
    '<rect x="132" y="284" width="52" height="32" fill="#2E7D46" opacity="0.95"/>'
    '<rect x="186" y="284" width="52" height="32" fill="#4C9A63" opacity="0.75"/>'
    '<rect x="240" y="284" width="52" height="32" fill="#8FCB9F" opacity="0.6"/>'
    '<rect x="294" y="284" width="52" height="32" fill="#C9E3D1" opacity="0.5"/>'
    '<rect x="348" y="284" width="52" height="32" fill="#C9E3D1" opacity="0.45"/>'
    '<rect x="402" y="284" width="54" height="32" fill="#8FCB9F" opacity="0.6"/>'
    '</g>'
    '</svg>'
)


def render_hero_illustration() -> None:
    """Render a stylised, non-data illustrative preview for the home hero."""

    st.markdown(_HERO_ILLUSTRATION_SVG, unsafe_allow_html=True)


def render_chip_row(labels: list[str]) -> None:
    """Render a wrapped row of pill-style capability chips."""

    chips = "".join(f'<span class="pee-chip">{label}</span>' for label in labels)
    st.html(f'<div class="pee-chip-row">{chips}</div>')
