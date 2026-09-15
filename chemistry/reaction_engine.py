"""
OrgReact — Reaction Profile Selector
"""

from dataclasses import dataclass

from chemistry.classification import classify_molecule
from chemistry.reagents import find_reagent
from chemistry.reaction_profiles import find_profiles, ReactionProfile


@dataclass(frozen=True)
class ReactionMatch:
    profile: ReactionProfile
    confidence: str
    reason: str


def get_molecule_labels(molecule) -> set[str]:
    classification = classify_molecule(molecule)

    labels = set()

    # Carbon labels
    labels.update(classification.carbon_types.values())

    # Alcohol labels
    for alcohol in classification.alcohol_sites:
        labels.add(alcohol["classification"])

    # Functional-group labels
    for group in classification.functional_groups:
        labels.add(group.name)

    return labels


def reagent_matches(profile: ReactionProfile, reagent_name: str) -> bool:
    """
    Handles profiles such as:

        KMnO4
        KMnO4/K2Cr2O7
        PCC/KMnO4/K2Cr2O7
        HCl/HBr/HI
    """

    requested = reagent_name.lower().replace(" ", "")
    profile_reagent = profile.reagent.lower().replace(" ", "")

    # Direct complete match
    if requested == profile_reagent:
        return True

    # Profile contains several alternative reagents
    alternatives = profile_reagent.split("/")

    return requested in alternatives


def substrate_matches(
    profile: ReactionProfile,
    molecule_labels: set[str],
) -> tuple[bool, str]:

    # Exact functional-group match
    for substrate in profile.substrate:
        if substrate in molecule_labels:
            return True, f"Matched substrate: {substrate}"

    # Generic alcohol profile
    if "alcohol" in profile.substrate:
        alcohol_types = {
            "primary alcohol",
            "secondary alcohol",
            "tertiary alcohol",
        }

        if molecule_labels.intersection(alcohol_types):
            return True, "Matched substrate: alcohol"

    return False, "No matching substrate"


def select_reactions(molecule, reagent_name: str):
    """
    Select reaction profiles using:

        molecule structure
        detected substrate
        reagent
    """

    reagent = find_reagent(reagent_name)

    if reagent is None:
        return {
            "success": False,
            "error": f"Unknown reagent: {reagent_name}",
            "labels": set(),
            "matches": [],
        }

    molecule_labels = get_molecule_labels(molecule)

    # Search all profiles first.
    candidate_profiles = find_profiles()

    matches = []

    for profile in candidate_profiles:

        if not reagent_matches(profile, reagent_name):
            continue

        substrate_ok, reason = substrate_matches(
            profile,
            molecule_labels,
        )

        if not substrate_ok:
            continue

        matches.append(
            ReactionMatch(
                profile=profile,
                confidence="high",
                reason=reason,
            )
        )

    return {
        "success": True,
        "reagent": reagent,
        "labels": molecule_labels,
        "matches": matches,
    }


def reaction_summary(result) -> str:
    if not result["success"]:
        return result["error"]

    lines = [
        "REACTION SELECTION",
        "=" * 40,
        f"Reagent: {result['reagent'].name}",
        "",
        "Detected molecular features:",
    ]

    for label in sorted(result["labels"]):
        lines.append(f"  - {label}")

    lines.append("")

    if not result["matches"]:
        lines.append("No matching reaction profile found.")
        return "\n".join(lines)

    lines.append("Matching reactions:")

    for match in result["matches"]:
        profile = match.profile

        lines.append("")
        lines.append(f"  {profile.id}")
        lines.append(f"  Reaction: {profile.reaction_type}")
        lines.append(f"  Conditions: {', '.join(profile.conditions)}")
        lines.append(f"  Product rule: {profile.product_rule}")
        lines.append(f"  Confidence: {match.confidence}")
        lines.append(f"  Reason: {match.reason}")

        if profile.notes:
            lines.append(f"  Notes: {profile.notes}")

    return "\n".join(lines)