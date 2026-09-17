"""Helpers for serializing VCF INFO fields."""


def format_vcf_info_item(key, value):
    """Serialize one INFO entry without corrupting VCF flag fields.

    Parsed INFO flags are represented internally as the boolean value ``True``.
    VCF flags must be emitted as the bare key (for example ``PRECISE``), not as
    ``PRECISE=True``.  Other values, including the literal string ``"True"``,
    remain ordinary key/value entries.
    """
    if value is True:
        return str(key)

    return f"{key}={value}"
