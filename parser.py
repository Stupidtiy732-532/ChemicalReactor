from molecule import Molecule


class Lexer:

    ELEMENTS = sorted(
        Molecule.ELEMENTS,
        key=len,
        reverse=True
    )

    def __init__(self, text):

        self.text = text.replace(" ", "")
        self.position = 0


    # ==================================================
    # TOKENIZER
    # ==================================================

    def tokenize(self):

        tokens = []

        while self.position < len(self.text):

            char = self.text[self.position]

            # ------------------------------------------
            # Parentheses
            # ------------------------------------------

            if char == "(":

                tokens.append(
                    ("LPAREN", "(")
                )

                self.position += 1
                continue


            if char == ")":

                tokens.append(
                    ("RPAREN", ")")
                )

                self.position += 1
                continue


            # ------------------------------------------
            # Bracket atom
            # ------------------------------------------

            if char == "[":

                tokens.append(
                    self.read_bracket_atom()
                )

                continue


            # ------------------------------------------
            # Radical
            # ------------------------------------------

            if char == ".":

                tokens.append(
                    ("RADICAL", ".")
                )

                self.position += 1
                continue


            # ------------------------------------------
            # Lone pair
            # ------------------------------------------

            if char == ":":

                tokens.append(
                    ("LONE_PAIR", ":")
                )

                self.position += 1
                continue


            # ------------------------------------------
            # Ring marker
            # ------------------------------------------

            if char == "@":

                self.position += 1

                number = self.read_number()

                if not number:

                    raise ValueError(
                        "Ring marker '@' must "
                        "be followed by a number."
                    )

                tokens.append(
                    ("RING", int(number))
                )

                continue


            # ------------------------------------------
            # Triple shorthand
            # -=
            # =-
            # ------------------------------------------

            if (
                self.position + 1
                < len(self.text)
            ):

                pair = (
                    self.text[
                        self.position:
                        self.position + 2
                    ]
                )

                if pair in ("-=", "=-"):

                    tokens.append(
                        ("BOND", 3)
                    )

                    self.position += 2
                    continue


            # ------------------------------------------
            # Bond symbols
            # ------------------------------------------

            if char == "-":

                tokens.append(
                    ("BOND", 1)
                )

                self.position += 1
                continue


            if char == "=":

                tokens.append(
                    ("BOND", 2)
                )

                self.position += 1
                continue


            if char in ("#", "≡"):

                tokens.append(
                    ("BOND", 3)
                )

                self.position += 1
                continue

            # ------------------------------------------
            # Ordinary atom
            # ------------------------------------------

            element = self.read_element()

            if element is not None:

                mass_number = None
                hydrogens = 0
                charge = 0
                marker = None

                # Isotope: C_13
                if (
                        self.position < len(self.text)
                        and self.text[self.position] == "_"
                ):

                    self.position += 1

                    number = self.read_number()

                    if not number:
                        raise ValueError(
                            "Mass marker '_' must "
                            "be followed by a number."
                        )

                    mass_number = int(number)

                # Hydrogens: H, H2, H3
                hydrogens = self.read_optional_hydrogens()

                # Charge: ^+, ^-, ^2+, ^+2
                charge = self.read_optional_charge()

                # Marker: <central> or X<central>
                if (
                        self.position < len(self.text)
                        and self.text[self.position] == "X"
                ):

                    self.position += 1

                    if (
                            self.position >= len(self.text)
                            or self.text[self.position] != "<"
                    ):
                        raise ValueError(
                            "Marker atom 'X' must be "
                            "followed by '<...>'."
                        )

                    marker_token = self.read_marker()
                    marker = marker_token[1]

                elif (
                        self.position < len(self.text)
                        and self.text[self.position] == "<"
                ):

                    marker_token = self.read_marker()
                    marker = marker_token[1]

                tokens.append(
                    ("ATOM", {
                        "element": element,
                        "hydrogens": hydrogens,
                        "charge": charge,
                        "radicals": 0,
                        "lone_pairs": 0,
                        "mass_number": mass_number,
                        "marker": marker
                    })
                )

                continue


            raise ValueError(
                f"Unexpected character "
                f"'{char}' at position "
                f"{self.position}."
            )


        tokens.append(
            ("EOF", None)
        )

        return tokens


    # ==================================================
    # ELEMENT
    # ==================================================

    def read_element(self):

        for element in self.ELEMENTS:

            if self.text.startswith(
                element,
                self.position
            ):

                self.position += len(element)

                return element

        return None


    # ==================================================
    # HYDROGENS
    # ==================================================

    def read_optional_hydrogens(self):

        if (
            self.position >= len(self.text)
            or self.text[self.position] != "H"
        ):

            return 0


        self.position += 1

        number = self.read_number()

        if number:

            return int(number)

        return 1

    # ==================================================
    # CHARGE
    #
    # ^+
    # ^-
    # ^2+
    # ^+2
    # ^++
    # ^--
    # ==================================================

    def read_optional_charge(self):

        if (
                self.position >= len(self.text)
                or self.text[self.position] != "^"
        ):
            return 0

        self.position += 1

        start = self.position

        # Digits before sign: ^2+
        leading_number = self.read_number()

        if (
                self.position >= len(self.text)
                or self.text[self.position] not in ("+", "-")
        ):
            raise ValueError(
                "Charge marker '^' must "
                "contain '+' or '-'."
            )

        sign = self.text[self.position]
        self.position += 1

        # Digits after sign: ^+2
        trailing_number = self.read_number()

        # Repeated signs: ^++, ^--
        repeated_signs = 1

        while (
                self.position < len(self.text)
                and self.text[self.position] == sign
        ):
            repeated_signs += 1
            self.position += 1

        number = (
                leading_number
                or trailing_number
        )

        magnitude = (
            int(number)
            if number
            else repeated_signs
        )

        if sign == "+":
            return magnitude

        return -magnitude

    # ==================================================
    # NUMBER
    # ==================================================

    def read_number(self):

        start = self.position

        while (
            self.position < len(self.text)
            and self.text[
                self.position
            ].isdigit()
        ):

            self.position += 1


        return (
            self.text[
                start:self.position
            ]
        )


    # ==================================================
    # MARKER
    #
    # <X>
    # <alpha>
    # <attack_site>
    # ==================================================

    def read_marker(self):

        self.position += 1

        start = self.position

        while (
            self.position < len(self.text)
            and self.text[
                self.position
            ] != ">"
        ):

            self.position += 1


        if self.position >= len(self.text):

            raise ValueError(
                "Unclosed atom marker '<...>'."
            )


        marker = self.text[
            start:self.position
        ]

        self.position += 1


        if not marker:

            raise ValueError(
                "Atom marker cannot be empty."
            )


        return (
            "MARKER",
            marker
        )


    # ==================================================
    # BRACKET ATOM
    # ==================================================

    def read_bracket_atom(self):

        self.position += 1

        start = self.position

        while (
            self.position < len(self.text)
            and self.text[
                self.position
            ] != "]"
        ):

            self.position += 1


        if self.position >= len(self.text):

            raise ValueError(
                "Unclosed bracket atom."
            )


        content = self.text[
            start:self.position
        ]

        self.position += 1


        if not content:

            raise ValueError(
                "Empty bracket atom."
            )


        return (
            "ATOM",
            self.parse_bracket_content(
                content
            )
        )


    # ==================================================
    # PARSE BRACKET CONTENT
    # ==================================================

    def parse_bracket_content(
        self,
        content
    ):

        position = 0

        element = None

        for candidate in self.ELEMENTS:

            if content.startswith(
                candidate,
                position
            ):

                element = candidate

                position += len(candidate)

                break


        if element is None:

            raise ValueError(
                f"Invalid bracket atom: "
                f"[{content}]"
            )


        hydrogens = 0
        charge = 0
        radicals = 0
        lone_pairs = 0
        mass_number = None
        marker = None


        while position < len(content):

            char = content[position]


            # Hydrogen

            if char == "H":

                position += 1

                start = position

                while (
                    position < len(content)
                    and content[position].isdigit()
                ):

                    position += 1


                number = content[
                    start:position
                ]

                hydrogens += (
                    int(number)
                    if number
                    else 1
                )

                continue


            # Mass number

            if char == "_":

                position += 1

                start = position

                while (
                    position < len(content)
                    and content[position].isdigit()
                ):

                    position += 1


                number = content[
                    start:position
                ]

                if not number:

                    raise ValueError(
                        "Mass marker '_' must "
                        "be followed by a number."
                    )


                mass_number = int(number)

                continue

            # Marker prefix: X<central>

            if content.startswith("X", position):

                position += 1

                if (
                        position >= len(content)
                        or content[position] != "<"
                ):
                    raise ValueError(
                        "Marker atom 'X' must be "
                        "followed by '<...>'."
                    )

                position += 1

                start = position

                while (
                        position < len(content)
                        and content[position] != ">"
                ):
                    position += 1

                if position >= len(content):
                    raise ValueError(
                        "Unclosed atom marker "
                        "inside bracket."
                    )

                marker = content[start:position]

                position += 1

                if not marker:
                    raise ValueError(
                        "Atom marker cannot be empty."
                    )

                continue


            # Marker

            if char == "<":

                position += 1

                start = position

                while (
                    position < len(content)
                    and content[position] != ">"
                ):

                    position += 1


                if position >= len(content):

                    raise ValueError(
                        "Unclosed atom marker "
                        "inside bracket."
                    )


                marker = content[
                    start:position
                ]

                position += 1

                if not marker:

                    raise ValueError(
                        "Atom marker cannot be empty."
                    )

                continue


            # Radical

            if char == ".":

                radicals += 1

                position += 1

                continue


            # Lone pair

            if char == ":":

                lone_pairs += 1

                position += 1

                continue


            # Charge notation
            #
            # ^+
            # ^-
            # ^2+
            # ^+2
            #

            if char == "^":

                position += 1

                start = position

                while (
                    position < len(content)
                    and content[position].isdigit()
                ):

                    position += 1


                leading_number = content[
                    start:position
                ]


                if (
                    position >= len(content)
                    or content[position]
                    not in ("+", "-")
                ):

                    raise ValueError(
                        "Charge marker '^' must "
                        "contain + or -."
                    )


                sign = content[position]

                position += 1


                start = position

                while (
                    position < len(content)
                    and content[position].isdigit()
                ):

                    position += 1


                trailing_number = content[
                    start:position
                ]


                number = (
                    leading_number
                    or trailing_number
                )


                magnitude = (
                    int(number)
                    if number
                    else 1
                )


                if sign == "+":

                    charge += magnitude

                else:

                    charge -= magnitude


                continue


            raise ValueError(
                f"Unknown modifier "
                f"'{char}' in "
                f"[{content}]"
            )


        return {
            "element": element,
            "hydrogens": hydrogens,
            "charge": charge,
            "radicals": radicals,
            "lone_pairs": lone_pairs,
            "mass_number": mass_number,
            "marker": marker
        }


# ==================================================
# PARSER
# ==================================================

class Parser:

    def __init__(self, text):

        self.tokens = Lexer(
            text
        ).tokenize()

        self.position = 0

        self.molecule = Molecule()

        self.current_atom = None

        # None means:
        # no explicit bond is pending

        self.pending_bond = None

        self.ring_connections = {}


    # ==================================================
    # TOKEN HELPERS
    # ==================================================

    def current(self):

        return self.tokens[
            self.position
        ]


    def advance(self):

        token = self.current()

        self.position += 1

        return token


    # ==================================================
    # MAIN PARSER
    # ==================================================

    def parse(self):

        self.parse_sequence()


        if self.current()[0] != "EOF":

            raise ValueError(
                "Unexpected input."
            )


        # Explicit dangling bond

        if self.pending_bond is not None:

            raise ValueError(
                "Structure ends with an incomplete bond."
            )


        if self.ring_connections:

            raise ValueError(
                "Unclosed ring connection(s): "
                + ", ".join(
                    str(x)
                    for x in self.ring_connections
                )
            )


        return self.molecule


    # ==================================================
    # SEQUENCE
    # ==================================================

    def parse_sequence(self):

        while True:

            token_type, value = self.current()


            if token_type in (
                "EOF",
                "RPAREN"
            ):

                return


            # Branch

            if token_type == "LPAREN":

                self.parse_branch()

                continue


            # Bond

            if token_type == "BOND":

                if self.pending_bond is not None:

                    raise ValueError(
                        "Two consecutive bond "
                        "operators."
                    )


                self.pending_bond = value

                self.advance()

                continue


            # Ring

            if token_type == "RING":

                self.handle_ring(
                    value
                )

                self.advance()

                continue


            # Prefix radical

            if token_type == "RADICAL":

                radical_count = 0


                while (
                    self.current()[0]
                    == "RADICAL"
                ):

                    radical_count += 1

                    self.advance()


                if self.current()[0] != "ATOM":

                    raise ValueError(
                        "Radical marker must "
                        "precede an atom."
                    )


                atom_data = dict(
                    self.current()[1]
                )

                atom_data[
                    "radicals"
                ] += radical_count


                self.parse_atom(
                    atom_data
                )

                continue


            # Atom

            if token_type == "ATOM":

                self.parse_atom(
                    dict(value)
                )

                continue


            raise ValueError(
                f"Unexpected token "
                f"{token_type}."
            )


    # ==================================================
    # ATOM
    # ==================================================

    def parse_atom(self, data):

        self.advance()


        # ----------------------------------------------
        # Post-atom metadata
        #
        # C_14
        # C<X>
        # C_14<X>
        # ----------------------------------------------

        while True:

            token_type, value = self.current()


            if token_type == "MASS":

                if data["mass_number"] is not None:

                    raise ValueError(
                        "Atom already has a "
                        "mass number."
                    )


                data["mass_number"] = value

                self.advance()

                continue


            if token_type == "MARKER":

                if data["marker"] is not None:

                    raise ValueError(
                        "Atom already has a "
                        "marker."
                    )


                data["marker"] = value

                self.advance()

                continue


            break


        atom_id = self.molecule.add_atom(

            element=data["element"],

            hydrogens=data[
                "hydrogens"
            ],

            charge=data[
                "charge"
            ],

            radical_electrons=data[
                "radicals"
            ],

            lone_pairs=data[
                "lone_pairs"
            ],

            mass_number=data[
                "mass_number"
            ],

            mark=data[
                "marker"
            ]
        )


        # Connect previous atom

        if self.current_atom is not None:

            bond_order = (
                self.pending_bond
                if self.pending_bond
                is not None
                else 1
            )


            self.molecule.add_bond(

                self.current_atom,

                atom_id,

                bond_order
            )


        self.current_atom = atom_id

        self.pending_bond = None


        # ----------------------------------------------
        # Postfix modifiers
        # ----------------------------------------------

        while True:

            token_type, value = self.current()


            # Radical

            if token_type == "RADICAL":

                atom = self.molecule.get_atom(
                    atom_id
                )

                atom.radical_electrons += 1

                self.advance()

                continue


            # Lone pair

            if token_type == "LONE_PAIR":

                atom = self.molecule.get_atom(
                    atom_id
                )

                atom.lone_pairs += 1

                self.advance()

                continue


            # Ring

            if token_type == "RING":

                self.handle_ring(
                    value
                )

                self.advance()

                continue


            break


    # ==================================================
    # BRANCH
    # ==================================================

    def parse_branch(self):

        self.advance()


        if self.current_atom is None:

            raise ValueError(
                "Branch cannot occur "
                "before an atom."
            )


        branch_origin = (
            self.current_atom
        )


        branch_bond = (
            self.pending_bond
        )


        self.pending_bond = None


        self.parse_sequence()


        if self.current()[0] != "RPAREN":

            raise ValueError(
                "Unclosed branch."
            )


        self.advance()


        self.current_atom = (
            branch_origin
        )


        self.pending_bond = (
            branch_bond
        )


    # ==================================================
    # RING
    # ==================================================

    def handle_ring(
        self,
        ring_number
    ):

        if self.current_atom is None:

            raise ValueError(
                "Ring marker must follow "
                "an atom."
            )


        if (
            ring_number
            not in self.ring_connections
        ):

            self.ring_connections[
                ring_number
            ] = {

                "atom": self.current_atom,

                "bond": self.pending_bond
            }


        else:

            connection = (
                self.ring_connections.pop(
                    ring_number
                )
            )


            other_atom = (
                connection["atom"]
            )


            stored_bond = (
                connection["bond"]
            )


            bond_order = (

                self.pending_bond

                if self.pending_bond
                is not None

                else (

                    stored_bond

                    if stored_bond
                    is not None

                    else 1
                )
            )


            self.molecule.add_bond(

                other_atom,

                self.current_atom,

                bond_order
            )


        self.pending_bond = None


# ==================================================
# PUBLIC FUNCTION
# ==================================================

def parse_structure(structure):

    parser = Parser(structure)

    return parser.parse()