# -*- coding: utf-8 -*-
"""
fence_scan — SHARED helper for the anti-fishing fence (a named project directive, S891).
================================================================================
★WHY THIS EXISTS (structural fix, memory-family item 12):
  the anti-fishing fence looks for a FORBIDDEN PHRASE (b₀=9-IN, «(3,1) обрана», MATCH-мінт)
  in a probe's text. But the very TEXT that discusses the phrase contains it ⟹ the fence
  catches its OWN author (S889: a visa quoted the forbidden phrase; the fixed regex
  caught its OWN regex line — recursively). Every probe rewrote the GUARDLINE
  filter by hand and COULD FORGET IT — Beta forgot it twice in one session.
  ⟹ the fix is NOT "remember to add GUARDLINE", but ONE place where excluding
  a discussion is IMPOSSIBLE to miss. Here it is BY CONSTRUCTION.

★INVARIANT: a line carrying the `# GUARDLINE` tag (or one of the extra `guard_markers`)
  is ALWAYS excluded from detection. A probe cannot "forget" — the exclusion is made
  by the helper, not the author.

Usage in a probe:
    from tools.fence_scan import scan_forbidden
    hits = scan_forbidden(__file__, [r"мінт\\w*\\s+b₀=9", r"обрано\\s+\\(3,1\\)"])
    assert not hits, f"fence: forbidden phrase outside GUARDLINE: {hits}"
"""
import re
from pathlib import Path

GUARD_MARK = "GUARDLINE"


def _clean_lines(text: str, guard_markers=()) -> str:
    """Strips discussion lines (GUARDLINE + extra markers) BY CONSTRUCTION."""
    marks = (GUARD_MARK,) + tuple(guard_markers)
    return "\n".join(ln for ln in text.split("\n")
                     if not any(m in ln for m in marks))


def scan_forbidden(source, patterns, guard_markers=(), ignorecase=True):
    """Returns the list of (pattern, matched_text) forbidden phrases found in
    `source`'s text OUTSIDE GUARDLINE. `source` = a path to a file OR the text
    itself. An empty list = the fence is clean (no mint outside a discussion).

    ★CASE INVARIANT (directive-2 per an internal project directive, from the W28/S908 visa
      finding, 2026-07-16): detection is CASE-INSENSITIVE BY DEFAULT.
      WHY: before the fix this was a bare re.finditer(pat, clean), and every
      probe's patterns were written in LOWERCASE — while the project marks its
      rulings in UPPERCASE ("P10 CONFIRMED BLIND — THE DIAGONAL IS DEAD",
      commit d53f71ee). ⟹ the fence covered half the surface, and the hole
      landed EXACTLY on the working case of rulings: "TIME = a ladder" ·
      "The pocket is STABLE" — both passed straight through.
      The blindness was held up by DISCIPLINE OF THE HAND, not by the fence
      (measured by mutation, S908).
      ⟹ invariant: a hand cannot dodge the fence with UPPERCASE letters. By
      construction. ignorecase=False is left only for deliberately
      case-sensitive targets (e.g. a symbol's name); it is NOT a silent
      default.
    """
    p = Path(source)
    text = p.read_text(encoding="utf-8") if p.exists() else str(source)
    clean = _clean_lines(text, guard_markers)
    flags = re.IGNORECASE if ignorecase else 0
    hits = []
    for pat in patterns:
        for m in re.finditer(pat, clean, flags):
            hits.append((pat, m.group(0)))
    return hits


def fence_ok(source, patterns, guard_markers=(), ignorecase=True):
    """True iff no forbidden phrase outside GUARDLINE. Convenient boolean form."""
    return not scan_forbidden(source, patterns, guard_markers, ignorecase)


# ── ★KILL-TEST (run as a script): the helper CAN fail on a real mint
#    and lets a discussion pass. Without this — "a test that cannot fail".
if __name__ == "__main__":
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    PAT = [r"мінт\w*\s+b₀=9", r"обрано\s+\(3,1\)"]
    cases = [
        # (text, expects a mint?)
        ("нейтральний рядок без мішені", False),
        ("тут я мінтую b₀=9 IN — справжній мінт", True),
        ("# обговорюю: не мінтую b₀=9  GUARDLINE", False),   # discussion under the tag
        ("рядок: обрано (3,1) як корінь", True),
        ("чому не можна обрано (3,1)  # GUARDLINE", False),  # discussion
        # ★CASE-KILL (directive-2; these lines used to pass straight through
        #   BEFORE the fix — exactly where the fence was blind, since rulings
        #   are marked in UPPERCASE):
        ("вирок: ОБРАНО (3,1) як корінь", True),            # uppercase GUARDLINE
        ("вирок: Обрано (3,1) як корінь", True),            # title case GUARDLINE
        ("тут я МІНТУЮ b₀=9 IN", True),                     # uppercase + \w*
        ("# чому не можна ОБРАНО (3,1)  GUARDLINE", False),  # uppercase UNDER the tag
    ]
    ok = True
    print("★KILL-TEST fence_scan (can fail on a mint, lets a discussion pass,")
    print("  and is NOT dodged by case — case invariant, directive-2/S909):")
    for text, expect_mint in cases:
        got = bool(scan_forbidden(text, PAT))
        good = got == expect_mint
        ok &= good
        print(f"  [{'✓' if good else '✗ FAIL'}] mint={got} (expected {expect_mint}) "
              f"· «{text[:46]}»")
    # recursive case: the pattern-definition LINE ITSELF carries the phrase, but is tagged
    recursive = 'PAT = [r"мінт b₀=9|обрано (3,1)"]  # GUARDLINE'
    rec_ok = not scan_forbidden(recursive, PAT)
    print(f"  [{'✓' if rec_ok else '✗'}] recursive: the definition line with the tag "
          f"does NOT catch itself — {rec_ok}")
    ok &= rec_ok
    print("ALL PASS" if ok else "★FAIL")
    sys.exit(0 if ok else 1)
