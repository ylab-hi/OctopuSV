"""
Quality filtering module for structural variants.

This module provides filtering capabilities based on variant-level attributes
such as quality scores, supporting reads, depth, and genotype quality.
"""

import logging
import re
from typing import Optional, Dict, Any


class QualityFilter:
    """Filter SVEvents based on quality attributes."""

    def __init__(self,
                 min_qual: Optional[float] = None,
                 max_qual: Optional[float] = None,
                 min_support: Optional[int] = None,
                 max_support: Optional[int] = None,
                 min_depth: Optional[int] = None,
                 max_depth: Optional[int] = None,
                 min_gq: Optional[int] = None,
                 min_svlen: Optional[int] = None,
                 max_svlen: Optional[int] = None,
                 filter_pass: bool = False,
                 exclude_nocall: bool = False):
        """Initialize quality filter with thresholds.

        Args:
            min_qual: Minimum QUAL score
            max_qual: Maximum QUAL score
            min_support: Minimum supporting reads
            max_support: Maximum supporting reads
            min_depth: Minimum total depth
            max_depth: Maximum total depth
            min_gq: Minimum genotype quality
            min_svlen: Minimum SV length
            max_svlen: Maximum SV length
            filter_pass: Only keep FILTER=PASS variants
            exclude_nocall: Exclude ./. genotypes
        """
        self.min_qual = min_qual
        self.max_qual = max_qual
        self.min_support = min_support
        self.max_support = max_support
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.min_gq = min_gq
        self.min_svlen = min_svlen
        self.max_svlen = max_svlen
        self.filter_pass = filter_pass
        self.exclude_nocall = exclude_nocall

        # Statistics
        self.stats = {
            'total': 0,
            'passed': 0,
            'filtered_qual': 0,
            'filtered_support': 0,
            'filtered_depth': 0,
            'filtered_gq': 0,
            'filtered_svlen': 0,
            'filtered_pass': 0,
            'filtered_nocall': 0
        }

    def filter_event(self, event) -> bool:
        """Check if an SVEvent passes all filter criteria.

        Args:
            event: SVEvent object

        Returns:
            True if event passes all filters, False otherwise
        """
        self.stats['total'] += 1

        # QUAL filter
        if not self._check_qual_filter(event):
            self.stats['filtered_qual'] += 1
            return False

        # FILTER field filter
        if not self._check_filter_field(event):
            self.stats['filtered_pass'] += 1
            return False

        # Support reads filter
        if not self._check_support_filter(event):
            self.stats['filtered_support'] += 1
            return False

        # Depth filter
        if not self._check_depth_filter(event):
            self.stats['filtered_depth'] += 1
            return False

        # Genotype quality filter
        if not self._check_gq_filter(event):
            self.stats['filtered_gq'] += 1
            return False

        # SV length filter
        if not self._check_svlen_filter(event):
            self.stats['filtered_svlen'] += 1
            return False

        # No-call genotype filter
        if not self._check_nocall_filter(event):
            self.stats['filtered_nocall'] += 1
            return False

        self.stats['passed'] += 1
        return True

    def _check_qual_filter(self, event) -> bool:
        """Check QUAL score filters."""
        if self.min_qual is None and self.max_qual is None:
            return True

        qual_value = self._get_qual_value(event)
        if qual_value is None:
            return True  # Skip filtering if QUAL is missing

        if self.min_qual is not None and qual_value < self.min_qual:
            return False
        if self.max_qual is not None and qual_value > self.max_qual:
            return False

        return True

    def _check_filter_field(self, event) -> bool:
        """Check FILTER field."""
        if not self.filter_pass:
            return True

        return event.filter.upper() == "PASS"

    def _check_support_filter(self, event) -> bool:
        """Check supporting reads filters."""
        if self.min_support is None and self.max_support is None:
            return True

        support_value = self._get_support_value(event)
        if support_value is None:
            return True  # Skip filtering if support info is missing

        if self.min_support is not None and support_value < self.min_support:
            return False
        if self.max_support is not None and support_value > self.max_support:
            return False

        return True

    def _check_depth_filter(self, event) -> bool:
        """Check depth filters."""
        if self.min_depth is None and self.max_depth is None:
            return True

        depth_value = self._get_depth_value(event)
        if depth_value is None:
            return True  # Skip filtering if depth info is missing

        if self.min_depth is not None and depth_value < self.min_depth:
            return False
        if self.max_depth is not None and depth_value > self.max_depth:
            return False

        return True

    def _check_gq_filter(self, event) -> bool:
        """Check genotype quality filter."""
        if self.min_gq is None:
            return True

        gq_value = self._get_gq_value(event)
        if gq_value is None:
            return True  # Skip filtering if GQ is missing

        return gq_value >= self.min_gq

    def _check_svlen_filter(self, event) -> bool:
        """Check SV length filters."""
        if self.min_svlen is None and self.max_svlen is None:
            return True

        svlen_value = self._get_svlen_value(event)
        if svlen_value is None:
            return True  # Skip filtering if SVLEN is missing

        if self.min_svlen is not None and svlen_value < self.min_svlen:
            return False
        if self.max_svlen is not None and svlen_value > self.max_svlen:
            return False

        return True

    def _check_nocall_filter(self, event) -> bool:
        """Check no-call genotype filter."""
        if not self.exclude_nocall:
            return True

        gt_value = self._get_gt_value(event)
        if gt_value is None:
            return True

        # Check for no-call patterns: ./., .|., .
        nocall_patterns = ["./.", ".|.", "."]
        return gt_value not in nocall_patterns

    def _get_qual_value(self, event) -> Optional[float]:
        """Extract QUAL value from event."""
        try:
            if event.qual == "." or event.qual is None:
                return None
            return float(event.qual)
        except (ValueError, AttributeError):
            return None

    def _get_support_value(self, event) -> Optional[int]:
        """Extract supporting reads value from event."""
        # Try INFO field first (standardized SUPPORT)
        if "SUPPORT" in event.info and event.info["SUPPORT"] != ".":
            return self._safe_int(event.info["SUPPORT"])

        # Try other INFO fields
        info_fields = ["SUPPREAD", "RE", "DV"]
        for field in info_fields:
            if field in event.info and event.info[field] != ".":
                return self._safe_int(event.info[field])

        # Try FORMAT/SAMPLE fields
        sample_dict = self._parse_sample_fields(event)
        sample_fields = ["DV", "DR", "AD"]

        for field in sample_fields:
            if field in sample_dict:
                value = sample_dict[field]
                if field == "AD":
                    # AD format: ref,alt - use alt count
                    parts = value.split(",")
                    if len(parts) >= 2:
                        return self._safe_int(parts[1])
                else:
                    return self._safe_int(value)

        return None

    def _get_depth_value(self, event) -> Optional[int]:
        """Extract depth value from event."""
        sample_dict = self._parse_sample_fields(event)

        # Try DP field first
        if "DP" in sample_dict:
            return self._safe_int(sample_dict["DP"])

        # Try to calculate from AD (ref + alt)
        if "AD" in sample_dict:
            ad_value = sample_dict["AD"]
            parts = ad_value.split(",")
            if len(parts) >= 2:
                ref_count = self._safe_int(parts[0])
                alt_count = self._safe_int(parts[1])
                if ref_count is not None and alt_count is not None:
                    return ref_count + alt_count

        return None

    def _get_gq_value(self, event) -> Optional[int]:
        """Extract genotype quality value from event."""
        sample_dict = self._parse_sample_fields(event)
        if "GQ" in sample_dict:
            return self._safe_int(sample_dict["GQ"])
        return None

    def _get_svlen_value(self, event) -> Optional[int]:
        """Extract SV length value from event."""
        if "SVLEN" in event.info and event.info["SVLEN"] != ".":
            svlen = self._safe_int(event.info["SVLEN"])
            return abs(svlen) if svlen is not None else None
        return None

    def _get_gt_value(self, event) -> Optional[str]:
        """Extract genotype value from event."""
        sample_dict = self._parse_sample_fields(event)
        return sample_dict.get("GT")

    def _parse_sample_fields(self, event) -> Dict[str, str]:
        """Parse FORMAT and SAMPLE fields into dictionary."""
        try:
            format_parts = event.format.split(":")
            sample_parts = event.sample.split(":")
            return dict(zip(format_parts, sample_parts, strict=False))
        except (AttributeError, ValueError):
            return {}

    def _safe_int(self, value) -> Optional[int]:
        """Safely convert value to integer."""
        if value is None or value == "." or value == "":
            return None
        try:
            return int(float(value))  # Handle "3.0" format
        except (ValueError, TypeError):
            return None

    def get_stats(self) -> Dict[str, Any]:
        """Get filtering statistics."""
        stats = self.stats.copy()
        if stats['total'] > 0:
            stats['pass_rate'] = stats['passed'] / stats['total']
        else:
            stats['pass_rate'] = 0.0
        return stats

    def print_stats(self):
        """Print filtering statistics."""
        stats = self.get_stats()
        logging.info(f"Quality filtering statistics:")
        logging.info(f"  Total variants: {stats['total']}")
        logging.info(f"  Passed filters: {stats['passed']}")
        logging.info(f"  Pass rate: {stats['pass_rate']:.2%}")

        if stats['filtered_qual'] > 0:
            logging.info(f"  Filtered by QUAL: {stats['filtered_qual']}")
        if stats['filtered_support'] > 0:
            logging.info(f"  Filtered by support: {stats['filtered_support']}")
        if stats['filtered_depth'] > 0:
            logging.info(f"  Filtered by depth: {stats['filtered_depth']}")
        if stats['filtered_gq'] > 0:
            logging.info(f"  Filtered by GQ: {stats['filtered_gq']}")
        if stats['filtered_svlen'] > 0:
            logging.info(f"  Filtered by SVLEN: {stats['filtered_svlen']}")
        if stats['filtered_pass'] > 0:
            logging.info(f"  Filtered by FILTER: {stats['filtered_pass']}")
        if stats['filtered_nocall'] > 0:
            logging.info(f"  Filtered by no-call GT: {stats['filtered_nocall']}")